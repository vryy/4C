/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fluid_fsi.H

\brief Fluid field adapter for fixed-grid/ALE-fluid

<pre>
Maintainer: Raffaela Kruse
            kruse@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#include "ad_fld_fluid_fluid_fsi.H"

#include "../drt_fluid_xfluid/xfluidfluid.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidFluidFSI::FluidFluidFSI(
    Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<DRT::Discretization> embfluiddis,
    Teuchos::RCP<DRT::Discretization> bgfluiddis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output,
    bool isale,
    bool dirichletcond)
  : FluidFSI(
      fluid,
      embfluiddis,
      solver,
      params,
      output,
      isale,
      dirichletcond)
{
  // cast fluid to XFluidFluid
  xfluidfluid_ = Teuchos::rcp_dynamic_cast<FLD::XFluidFluid>(fluid_);
  if (xfluidfluid_ == Teuchos::null)
    dserror("Failed to cast ADAPTER::Fluid to FLD::XFluidFluid.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::Init()
{
  // call base class init
  FluidWrapper::Init();

  // determine the type of monolithic approach
  const Teuchos::ParameterList& xfluiddyn  = params_->sublist("XFLUID DYNAMIC/GENERAL");
  monolithic_approach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>
               (xfluiddyn,"MONOLITHIC_XFFSI_APPROACH");

  // should ALE-relaxation be carried out?
  relaxing_ale_ = (bool)DRT::INPUT::IntegralValue<int>(xfluiddyn,"RELAXING_ALE");

  // get no. of timesteps, after which ALE-mesh should be relaxed
  relaxing_ale_every_ = xfluiddyn.get<int>("RELAXING_ALE_EVERY");

  if (! relaxing_ale_ && relaxing_ale_every_ != 0)
    dserror("You don't want to relax the ALE but provide a relaxation interval != 0 ?!");

  if (relaxing_ale_every_ < 0)
    dserror("Please provide a reasonable relaxation interval. We can't travel back in time yet.");

  // create map extractor for combined fluid domains
  // (to distinguish between FSI interface DOF / merged inner embedded & background fluid DOF)
  mergedfluidinterface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  SetupInterface();

  interfaceforcen_ = Teuchos::rcp(new Epetra_Vector(*(FluidFSI::Interface()->FSICondMap())));

  // build inner velocity dof map (no DBC, no FSI)
  FluidFSI::BuildInnerVelMap();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::PrepareTimeStep()
{
  xfluidfluid_->PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFluidFSI::DofRowMap()
{
  return xfluidfluid_->DofRowMap();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::Solve()
{
  // cut and do XFEM time integration
  xfluidfluid_->PrepareSolve();
  // solve
  xfluidfluid_->Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::Update()
{
  if ( IsAleRelaxationStep(Step())
       && (monolithic_approach_ == INPAR::XFEM::XFFSI_FixedALE_Partitioned ||
           monolithic_approach_ == INPAR::XFEM::XFFSI_FixedALE_Interpolation))
  {
    // cut to abtain new interface position
    xfluidfluid_->CutAndSetState();

    // perform xfem time integration
    xfluidfluid_->DoTimeStepTransfer();

    if (monolithic_approach_ == INPAR::XFEM::XFFSI_FixedALE_Partitioned)
      xfluidfluid_->UpdateMonolithicFluidSolution(FluidFSI::Interface()->FSICondMap());

    // refresh the merged fluid map extractor
    SetupInterface();
    // create new extended shape derivatives matrix
    PrepareShapeDerivatives();
  }

  FluidWrapper::Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::XFluidFluidMapExtractor>const& ADAPTER::FluidFluidFSI::XFluidFluidMapExtractor()
{
  return xfluidfluid_->XFluidFluidMapExtractor();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::ApplyEmbFixedMeshDisplacement(Teuchos::RCP<const Epetra_Vector> disp)
{
  if (meshmap_ == Teuchos::null)
    dserror("Uninitialized mesh map");
  meshmap_->InsertCondVector(disp,xfluidfluid_->WriteAccessDispOldState());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  // it transfers the displacement we get from Ale-dis to the displacement of the
  // embedded-fluid-dis
  if (meshmap_ == Teuchos::null)
    dserror("Uninitialized mesh map");

  meshmap_->InsertCondVector(fluiddisp, xfluidfluid_->WriteAccessDispnp());

  // new grid velocity
  xfluidfluid_->UpdateGridv();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  if (meshmap_ == Teuchos::null)
    dserror("Missing mesh map!");
  meshmap_->InsertCondVector(gridvel,xfluidfluid_->WriteAccessGridVel());
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::FluidFluidFSI::BlockSystemMatrix()
{
  if (mergedfluidinterface_ == Teuchos::null)
    dserror("Uninitialized map FSI/inner fluid map extractor! Failed to create fluid block matrix.");

  // Create a local copy of the inner & conditioned map
  // Reason: the matrix splitting method from LINALG expects non-const maps
  Teuchos::RCP<Epetra_Map> innermap = Teuchos::rcp(new Epetra_Map(*mergedfluidinterface_->OtherMap()));
  Teuchos::RCP<Epetra_Map> condmap = Teuchos::rcp(new Epetra_Map(*mergedfluidinterface_->FSICondMap()));
  return xfluidfluid_->BlockSystemMatrix(innermap, condmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
)
{
  // call the usual routine
  xfluidfluid_->Evaluate(stepinc);

  // for fixed ALE approach, we only refresh the global fluid map extractor in Update()
  if (monolithic_approach_ != INPAR::XFEM::XFFSI_Full_Newton)
    return;

  // this is the case of a full Newton approach: update the map extractor, as fluid DOFs possibly have changed!
  SetupInterface();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFluidFSI::VelocityRowMap()
{
  // in case of fsi with fluidsplit, return the embedded velocity DOF
  // (to understand the motivation behind this, have a look at the recovery of the
  // Lagrange multiplier in standard ALE-FSI class (fluidsplit) in case of active
  // shape derivatives)
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(xfluidfluid_->XFluidFluidMapExtractor()->FluidMap());
  maps.push_back(xfluidfluid_->VelocityRowMap());
  Teuchos::RCP<const Epetra_Map> innervelocitymap = LINALG::MultiMapExtractor::IntersectMaps(maps);
  return innervelocitymap;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<DRT::Discretization>& ADAPTER::FluidFluidFSI::Discretization()
{
  return xfluidfluid_->EmbeddedDiscretization();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<IO::DiscretizationWriter>& ADAPTER::FluidFluidFSI::DiscWriter()
{
  return xfluidfluid_->EmbeddedDiscWriter();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::UseBlockMatrix(bool split_fluidsysmat)
{
  PrepareShapeDerivatives();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool ADAPTER::FluidFluidFSI::IsAleRelaxationStep(int step) const
{
  return relaxing_ale_ && step % relaxing_ale_every_ == 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::SetupInterface()
{
  if (mergedfluidinterface_ == Teuchos::null)
  {
    std::stringstream errmsg;
    errmsg
      << "Uninitialized map  map extractor for merged background & embedded inner/FSI fluid DOFs."
      << "\nFailed to perform map extractor setup.";
    dserror(errmsg.str());
  }

  FluidFSI::SetupInterface();

  // get background fluid map
  Teuchos::RCP<const Epetra_Map> xfluidmap = xfluidfluid_->XFluidFluidMapExtractor()->XFluidMap();
  // do the setup
  mergedfluidinterface_->Setup(xfluidmap, *FluidFSI::Interface());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::PrepareShapeDerivatives()
{
  // the dof-maps may have changed: create a new shape derivatives matrix
  Teuchos::RCP<std::set<int> > condelements = mergedfluidinterface_->ConditionedElementMap(*xfluidfluid_->EmbeddedDiscretization());
  xfluidfluid_->PrepareShapeDerivatives(mergedfluidinterface_,condelements);
}

