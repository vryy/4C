/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fluid_fsi.cpp

\brief Fluid field adapter for embedded (ALE-)fluid-fluid problems using XFEM

\level 2

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
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
ADAPTER::FluidFluidFSI::FluidFluidFSI(Teuchos::RCP<Fluid> xfluidfluid, Teuchos::RCP<Fluid> embfluid,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params, bool isale,
    bool dirichletcond)
    : FluidFSI(xfluidfluid, embfluid->Discretization(), solver, params,
          embfluid->Discretization()->Writer(), isale, dirichletcond)
{
  // cast fluid to XFluidFluid
  xfluidfluid_ = Teuchos::rcp_dynamic_cast<FLD::XFluidFluid>(xfluidfluid);
  if (xfluidfluid_ == Teuchos::null) dserror("Failed to cast ADAPTER::Fluid to FLD::XFluidFluid.");
  fluidimpl_ = Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(embfluid);
  if (fluidimpl_ == Teuchos::null)
    dserror("Failed to cast ADAPTER::Fluid to FLD::FluidImplicitTimInt.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::Init()
{
  // determine the type of monolithic approach
  const Teuchos::ParameterList& xfluiddyn = params_->sublist("XFLUID DYNAMIC/GENERAL");
  monolithic_approach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>(
      xfluiddyn, "MONOLITHIC_XFFSI_APPROACH");

  // should ALE-relaxation be carried out?
  relaxing_ale_ = (bool)DRT::INPUT::IntegralValue<int>(xfluiddyn, "RELAXING_ALE");

  // get no. of timesteps, after which ALE-mesh should be relaxed
  relaxing_ale_every_ = xfluiddyn.get<int>("RELAXING_ALE_EVERY");

  if (!relaxing_ale_ && relaxing_ale_every_ != 0)
    dserror("You don't want to relax the ALE but provide a relaxation interval != 0 ?!");

  if (relaxing_ale_every_ < 0)
    dserror("Please provide a reasonable relaxation interval. We can't travel back in time yet.");

  // create map extractor for combined fluid domains
  // (to distinguish between FSI interface DOF / merged inner embedded & background fluid DOF)
  mergedfluidinterface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  // call base class init
  FluidFSI::Init();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::PrepareTimeStep()
{
  if (Interface()->FSICondRelevant() &&
      (monolithic_approach_ == INPAR::XFEM::XFFSI_FixedALE_Partitioned ||
          monolithic_approach_ == INPAR::XFEM::XFFSI_FixedALE_Interpolation))
  {
    xfluidfluid_->SetInterfaceFixed();
  }
  else
  {
    xfluidfluid_->SetInterfaceFree();
  }
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
  // cut and do XFEM time integration, solve
  xfluidfluid_->Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::Update()
{
  if (Interface()->FSICondRelevant() && IsAleRelaxationStep(Step()) &&
      (monolithic_approach_ == INPAR::XFEM::XFFSI_FixedALE_Partitioned ||
          monolithic_approach_ == INPAR::XFEM::XFFSI_FixedALE_Interpolation))
  {
    // allow new interface position
    xfluidfluid_->SetInterfaceFree();

    // cut with new interface location and do XFEM time integration
    xfluidfluid_->PrepareXFEMSolve();

    // fix interface position again
    xfluidfluid_->SetInterfaceFixed();

    if (monolithic_approach_ == INPAR::XFEM::XFFSI_FixedALE_Partitioned)
      xfluidfluid_->UpdateMonolithicFluidSolution(FluidFSI::Interface()->FSICondMap());

    if (monolithic_approach_ == INPAR::XFEM::XFFSI_FixedALE_Interpolation)
      xfluidfluid_->InterpolateEmbeddedStateVectors();

    // refresh the merged fluid map extractor
    SetupInterface();

    // create new extended shape derivatives matrix
    PrepareShapeDerivatives();
  }

  FluidWrapper::Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::XFluidFluidMapExtractor> const&
ADAPTER::FluidFluidFSI::XFluidFluidMapExtractor()
{
  return xfluidfluid_->XFluidFluidMapExtractor();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  // store old state
  Teuchos::RCP<const Epetra_Vector> disp = meshmap_->ExtractCondVector(fluidimpl_->Dispnp());
  meshmap_->InsertCondVector(disp, xfluidfluid_->WriteAccessDispOldState());
  // apply mesh displacement and update grid velocity
  FluidFSI::ApplyMeshDisplacement(fluiddisp);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::FluidFluidFSI::BlockSystemMatrix()
{
  if (mergedfluidinterface_ == Teuchos::null)
    dserror(
        "Uninitialized map FSI/inner fluid map extractor! Failed to create fluid block matrix.");

  // Create a local copy of the inner & conditioned map
  // Reason: the matrix splitting method from LINALG expects non-const maps
  Teuchos::RCP<Epetra_Map> innermap =
      Teuchos::rcp(new Epetra_Map(*mergedfluidinterface_->OtherMap()));
  Teuchos::RCP<Epetra_Map> condmap =
      Teuchos::rcp(new Epetra_Map(*mergedfluidinterface_->FSICondMap()));
  return xfluidfluid_->BlockSystemMatrix(innermap, condmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::Evaluate(
    Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
)
{
  if (monolithic_approach_ == INPAR::XFEM::XFFSI_Full_Newton)
    *xfluidfluid_->WriteAccessDispOldState() = *fluidimpl_->Dispnp();

  // call the usual routine
  xfluidfluid_->Evaluate(stepinc);

  // for fixed ALE approach, we only refresh the global fluid map extractor in Update()
  if (monolithic_approach_ != INPAR::XFEM::XFFSI_Full_Newton) return;

  // this is the case of a full Newton approach: update the map extractor, as fluid DOFs possibly
  // have changed!
  SetupInterface();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidFSI::GridVel()
{
  return fluidimpl_->GridVel();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFluidFSI::WriteAccessGridVel()
{
  return fluidimpl_->WriteAccessGridVel();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidFSI::Dispnp() { return fluidimpl_->Dispnp(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFluidFSI::WriteAccessDispnp()
{
  return fluidimpl_->WriteAccessDispnp();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidFSI::Dispn() { return fluidimpl_->Dispn(); }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<DRT::Discretization>& ADAPTER::FluidFluidFSI::Discretization()
{
  return fluidimpl_->Discretization();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFluidFSI::VelocityRowMap()
{
  // in case of fsi with fluidsplit, return the embedded velocity DOF
  // (to understand the motivation behind this, have a look at the recovery of the
  // Lagrange multiplier in standard ALE-FSI class (fluidsplit) in case of active
  // shape derivatives)
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(xfluidfluid_->XFluidFluidMapExtractor()->FluidMap());
  maps.push_back(xfluidfluid_->VelocityRowMap());
  Teuchos::RCP<const Epetra_Map> innervelocitymap = LINALG::MultiMapExtractor::IntersectMaps(maps);
  return innervelocitymap;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::UseBlockMatrix(bool split_fluidsysmat) { PrepareShapeDerivatives(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool ADAPTER::FluidFluidFSI::IsAleRelaxationStep(int step) const
{
  return relaxing_ale_ && step % relaxing_ale_every_ == 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::FluidFluidFSI::ShapeDerivatives()
{
  return xfluidfluid_->ExtendedShapeDerivatives();
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
  Teuchos::RCP<std::set<int>> condelements =
      mergedfluidinterface_->ConditionedElementMap(*fluidimpl_->Discretization()());
  xfluidfluid_->PrepareShapeDerivatives(mergedfluidinterface_, condelements);
}
