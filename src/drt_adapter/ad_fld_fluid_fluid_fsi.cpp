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

#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_fluid_xfluid/xfluidfluid.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <vector>
#include <set>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidFluidFSI::FluidFluidFSI(Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<DRT::Discretization> embfluiddis,
    Teuchos::RCP<DRT::Discretization> bgfluiddis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    bool isale,
    bool dirichletcond,
    bool monolithicfluidfluidfsi)
  : FluidWrapper(fluid),
    embfluiddis_(embfluiddis),
    bgfluiddis_(bgfluiddis),
    solver_(solver),
    params_(params),
    dirichletcond_(dirichletcond),
    monolithicfluidfluidfsi_(monolithicfluidfluidfsi)
{
  // make sure
  if (fluid_ == Teuchos::null)
    dserror("Failed to create the underlying fluid adapter");

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::Init()
{
  // call base class init
  FluidWrapper::Init();

  // cast fluid to XFluidFluid
  xfluidfluid_ = Teuchos::rcp_dynamic_cast<FLD::XFluidFluid>(fluid_);
  if (xfluidfluid_ == Teuchos::null)
    dserror("Failed to cast ADAPTER::Fluid to FLD::XFluidFluid.");

  monolithic_approach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>
                        (params_->sublist("XFLUID DYNAMIC/GENERAL"),"MONOLITHIC_XFFSI_APPROACH");

  // map extractor for embedded fluid discretization
  interface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  interface_->Setup(*embfluiddis_);

  // map extractor for transfer of ALE-displacements to embedded discretization
  meshmap_ = Teuchos::rcp(new LINALG::MapExtractor());

  // map extractor for both fluid domains
  // (to distinguish between FSI interface DOF / merged inner embedded & background fluid DOF)
  mergedfluidinterface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  SetupInterfaceExtractor();

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint

  // here we get the dirichletmaps for the both discretizations
  const Teuchos::RCP<const LINALG::MapExtractor> embdbcmaps = xfluidfluid_->EmbeddedDirichMaps();
  const Teuchos::RCP<const LINALG::MapExtractor> bgdbcmaps = xfluidfluid_->BackgroundDirichMaps();

  // first build the inner map of embedded fluid (other map)
  // intersected with the dofs with no dbc
  {
    std::vector<Teuchos::RCP<const Epetra_Map> > maps;
    maps.push_back(interface_->OtherMap());
    maps.push_back(embdbcmaps->OtherMap());
    Teuchos::RCP<Epetra_Map> innervelmap_emb = LINALG::MultiMapExtractor::IntersectMaps(maps);

    // now the not-dbc map of background fluid and merge it with the
    // inner map of embedded fluid
    std::vector<Teuchos::RCP<const Epetra_Map> > bgembmaps;
    bgembmaps.push_back(bgdbcmaps->OtherMap());
    bgembmaps.push_back(innervelmap_emb);
    Teuchos::RCP<Epetra_Map> innermap_bgemb = LINALG::MultiMapExtractor::MergeMaps(bgembmaps);

    // now throw out the pressure dofs
    std::vector<Teuchos::RCP<const Epetra_Map> > finalmaps;
    finalmaps.push_back(innermap_bgemb);
    finalmaps.push_back(xfluidfluid_->VelocityRowMap());
    innervelmap_ = LINALG::MultiMapExtractor::IntersectMaps(finalmaps);
  }

  if (dirichletcond_)
  {
    // mark all interface velocities as dirichlet values
    xfluidfluid_->AddDirichCond(interface_->FSICondMap());
  }

  interfaceforcen_ = Teuchos::rcp(new Epetra_Vector(*(interface_->FSICondMap())));

  isfluidsplit_ = false;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidFluidFSI::TimeScaling() const
{
  if (params_->get<bool>("interface second order"))
    return 2./xfluidfluid_->Dt();
  else
    return 1./xfluidfluid_->Dt();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::Update()
{
  Teuchos::RCP<Epetra_Vector> interfaceforcem = interface_->ExtractFSICondVector(xfluidfluid_->TrueResidual());

  interfaceforcen_ = xfluidfluid_->ExtrapolateEndPoint(interfaceforcen_,interfaceforcem);

  xfluidfluid_->TimeUpdate();

  if (monolithicfluidfluidfsi_)
  {
    // refresh the merged fluid map extractor
    SetupInterfaceExtractor();
    // create new extended shape derivatives matrix
    PrepareShapeDerivatives();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFluidFSI::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFluidFSI::ExtractInterfaceForces()
{
  Teuchos::RCP<Epetra_Vector> interfaceforcem = interface_->ExtractFSICondVector(xfluidfluid_->TrueResidual());

  return xfluidfluid_->ExtrapolateEndPoint(interfaceforcen_,interfaceforcem);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFluidFSI::ExtractInterfaceVelnp()
{
  return interface_->ExtractFSICondVector(xfluidfluid_->Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFluidFSI::ExtractInterfaceVeln()
{
  return interface_->ExtractFSICondVector(xfluidfluid_->Veln());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor>const& ADAPTER::FluidFluidFSI::XFluidFluidMapExtractor()
{
  return xfluidfluid_->XFluidFluidMapExtractor();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::ApplyEmbFixedMeshDisplacement(Teuchos::RCP<const Epetra_Vector> disp)
{
  if (meshmap_ == Teuchos::null)
    dserror("Uninitialized mesh map");
  meshmap_->InsertCondVector(disp,xfluidfluid_->ViewOfDispoldstate());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  // it transfers the displacement we get from Ale-dis to the displacement of the
  // embedded-fluid-dis
  if (meshmap_ == Teuchos::null)
    dserror("Uninitialized mesh map");
  meshmap_->InsertCondVector(fluiddisp, xfluidfluid_->ViewOfDispnp());

  // new grid velocity
  xfluidfluid_->UpdateGridv();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  // apply the interface velocities
  interface_->InsertFSICondVector(ivel,xfluidfluid_->ViewofVelnp());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_->Setup(*embfluiddis_->DofRowMap(),mm,LINALG::SplitMap(*embfluiddis_->DofRowMap(),*mm));
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  if (meshmap_ == Teuchos::null)
    dserror("Missing mesh map!");
  meshmap_->InsertCondVector(gridvel,xfluidfluid_->ViewOfGridVel());
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = interface_->ExtractFSICondVector(Veln());
  /// We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  /*
   * Delta u(n+1,i+1) = fac * Delta d(n+1,i+1) - dt * u(n)
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  double timescale = TimeScaling();
  fcx->Update(-timescale*xfluidfluid_->Dt(),*veln,timescale);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = interface_->ExtractFSICondVector(Veln());

  /*
   * Delta d(n+1,i+1) = fac * [Delta u(n+1,i+1) + 2 * u(n)]
   *
   *             / = dt / 2   if interface time integration is second order
   * with fac = |
   *             \ = dt       if interface time integration is first order
   */
  double timescale = 1./TimeScaling();
  fcx->Update(xfluidfluid_->Dt(),*veln,timescale);
}

/*----------------------------------------------------------------------*
 * Remove passed DOFs from Dirichlet map (required for the monolithic
 * fluid-fluid fluidsplit algorithm)
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  xfluidfluid_->RemoveDirichCond(maptoremove);
  return;
}

/*----------------------------------------------------------------------*
 * Returns the embedded fluid DBC-MapExtractor() in the case of
 * monolithic fluid-split
 *----------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MapExtractor> ADAPTER::FluidFluidFSI::GetDBCMapExtractor()
{
  return xfluidfluid_->EmbeddedDirichMaps();
}

/*----------------------------------------------------------------------*
 *  method needed for block matrix creation on the fluid side!
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::FluidFluidFSI::BlockSystemMatrix()
{
  if (mergedfluidinterface_ == Teuchos::null)
    dserror("Uninitialized map FSI/inner fluid map extractor! Failed to create fluid block matrix.");

  // Create a local copy of the inner & conditioned map
  // Reason: the matrix splitting method from LINALG expects non-const maps
  Teuchos::RCP<Epetra_Map> innermap = Teuchos::rcp(new Epetra_Map(*mergedfluidinterface_->OtherMap()));
  Teuchos::RCP<Epetra_Map> condmap = Teuchos::rcp(new Epetra_Map(*mergedfluidinterface_->FSICondMap()));
  return xfluidfluid_->BlockSystemMatrix(mergedfluidinterface_, innermap, condmap);
}
/*---------------------------------------------------------------
* In case of fluid split:
* For correct output of the Lagrange Multiplier field fsilambda_,
* which belongs to the embedded fluid discretization
* in case of fluidsplit FSI, we need the
* DiscretizationWriter for the embedded fluid.
* (In default, DiscWriter() returns the DiscretizationWriter for
* the background fluid discretization!)
* We use the flag isfluidsplit_ for this.
* It is set to 'true' in UseBlockMatrix(bool), which is called with true
* in case of fluidsplit.
*------------------------------------------------------------------------*/
const Teuchos::RCP<IO::DiscretizationWriter>& ADAPTER::FluidFluidFSI::DiscWriter()
{
  if (isfluidsplit_)
    return xfluidfluid_->EmbDiscWriter();
  else
    return xfluidfluid_->DiscWriter();
}

void ADAPTER::FluidFluidFSI::Evaluate(
    Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
    )
{
  // call the usual routine
  xfluidfluid_->Evaluate(stepinc);

  // for fixed ALE approach, we only refresh the global fluid map extractor in Update()
  if (! monolithicfluidfluidfsi_ || monolithic_approach_ != INPAR::XFEM::XFFSI_Full_Newton)
    return;

  // this is the case of a full Newton approach: update the map extractor, as fluid DOFs possibly have changed!
  SetupInterfaceExtractor();
}

Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFluidFSI::VelocityRowMap()
{
  if (! monolithicfluidfluidfsi_ || ! isfluidsplit_)
    return fluid_->VelocityRowMap();

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

void ADAPTER::FluidFluidFSI::SetupInterfaceExtractor()
{
  if (mergedfluidinterface_ == Teuchos::null)
  {
    std::stringstream errmsg;
    errmsg
      << "Uninitialized map  map extractor for merged background & embedded inner/FSI fluid DOFs."
      << "\nFailed to perform map extractor setup.";
    dserror(errmsg.str());
  }

  if (interface_ == Teuchos::null)
  {
    std::stringstream errmsg;
    errmsg
      << "Uninitialized map extractor for solely embedded inner/FSI fluid DOFs."
      << "\nFailed to setup map extractor for merged background & embedded inner/FSI fluid DOFs.";
    dserror(errmsg.str());
  }

  // get background fluid map
  Teuchos::RCP<const Epetra_Map> xfluidmap = xfluidfluid_->XFluidFluidMapExtractor()->XFluidMap();
  // do the setup
  mergedfluidinterface_->Setup(xfluidmap, *interface_);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::PrepareShapeDerivatives()
{
  if (! monolithicfluidfluidfsi_)
    return;

  // the dof-maps may have changed: create a new shape derivatives matrix
  Teuchos::RCP<std::set<int> > condelements = interface_->ConditionedElementMap(*Discretization());
  xfluidfluid_->PrepareShapeDerivatives(mergedfluidinterface_,condelements);
}

/*----------------------------------------------------------------------*
 * request fluid system matrix as block matrix (fluidsplit FSI)
 * and allocate a new shape derivatives matrix
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFluidFSI::UseBlockMatrix(bool split_fluidsysmat)
{
  isfluidsplit_ = split_fluidsysmat;
  PrepareShapeDerivatives();
}
