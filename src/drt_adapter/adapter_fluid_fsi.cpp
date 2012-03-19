/*----------------------------------------------------------------------*/
#include "adapter_fluid_fsi.H"

#include "../drt_adapter/adapter_fluid.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../linalg/linalg_mapextractor.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <vector>
#include <set>
/*======================================================================*/
/* constructor */
ADAPTER::FluidFSI::FluidFSI(Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output,
    bool isale,
    bool dirichletcond)
: FluidWrapper(fluid),
  dis_(dis),
  solver_(solver),
  params_(params),
  output_(output)
{
  // make sure
  if (fluid_ == null)
    dserror("Failed to create the underlying fluid adapter");

  // cast fluid to fluidimplicit
  fluidimpl_ = Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_);
  if (fluidimpl_ == Teuchos::null)
    dserror("Failed to cast ADAPTER::Fluid to FLD::FluidImplicitTimeInt.");

  interface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  meshmap_   = Teuchos::rcp(new LINALG::MapExtractor());

  interface_->Setup(*dis);
  fluidimpl_->SetSurfaceSplitter(&(*interface_));

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint
  const Teuchos::RCP<const LINALG::MapExtractor> dbcmaps = fluidimpl_->DirichMaps();
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(interface_->OtherMap());
  maps.push_back(dbcmaps->OtherMap());
  innervelmap_ = LINALG::MultiMapExtractor::IntersectMaps(maps);

  if (dirichletcond)
  {
    // mark all interface velocities as dirichlet values
    fluidimpl_->AddDirichCond(interface_->FSICondMap());
  }

  interfaceforcen_ = rcp(new Epetra_Vector(*(interface_->FSICondMap())));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFSI::DofRowMap()
{
  return DofRowMap(0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFSI::DofRowMap(unsigned nds)
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap(nds);
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidFSI::TimeScaling() const
{
  if (params_->get<bool>("interface second order"))
    return 2./fluidimpl_->Dt();
  else
    return 1./fluidimpl_->Dt();
}


void ADAPTER::FluidFSI::Update()
{
  Teuchos::RCP<Epetra_Vector> interfaceforcem = interface_->ExtractFSICondVector(fluidimpl_->TrueResidual());

  interfaceforcen_ = fluidimpl_->ExtrapolateEndPoint(interfaceforcen_,interfaceforcem);

  fluidimpl_->TimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFSI::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractInterfaceForces()
{
  Teuchos::RCP<Epetra_Vector> interfaceforcem = interface_->ExtractFSICondVector(fluidimpl_->TrueResidual());

  return fluidimpl_->ExtrapolateEndPoint(interfaceforcen_,interfaceforcem);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractInterfaceFluidVelocity()
{
  return interface_->ExtractFSICondVector(fluidimpl_->Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractInterfaceVeln()
{
  return interface_->ExtractFSICondVector(fluidimpl_->Veln());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractFreeSurfaceVeln()
{
  return Interface().ExtractFSCondVector(fluidimpl_->Veln());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  meshmap_->InsertCondVector(fluiddisp,fluidimpl_->ViewOfDispnp());

  // new grid velocity
  fluidimpl_->UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  meshmap_->InsertCondVector(gridvel,fluidimpl_->ViewOfGridVel());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSICondVector(Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = TimeScaling();
  fcx->Update(-timescale*fluidimpl_->Dt(),*veln,timescale);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSICondVector(Veln());

  // We convert Delta u(n+1,i+1) to Delta d(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = 1./TimeScaling();
  fcx->Update(fluidimpl_->Dt(),*veln,timescale);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSCondVector(Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = TimeScaling();
  fcx->Update(-timescale*fluidimpl_->Dt(),*veln,timescale);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSCondVector(Veln());

  // We convert Delta u(n+1,i+1) to Delta d(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = 1./TimeScaling();
  fcx->Update(fluidimpl_->Dt(),*veln,timescale);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::IntegrateInterfaceShape()
{
  return interface_->ExtractFSICondVector(fluidimpl_->IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::UseBlockMatrix(bool splitmatrix)
{
  Teuchos::RCP<std::set<int> > condelements = Interface().ConditionedElementMap(*Discretization());
  fluidimpl_->UseBlockMatrix(condelements,Interface(),Interface(),splitmatrix);
}
