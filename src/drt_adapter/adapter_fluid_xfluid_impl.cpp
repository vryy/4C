/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_xfluid_impl.cpp


<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_fluid_xfluid_impl.H"
#include "../drt_lib/drt_condition_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidXFluidImpl::FluidXFluidImpl(
        Teuchos::RCP<DRT::Discretization> fluiddis,
        Teuchos::RCP<DRT::Discretization> xfluiddis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output,
        bool isale,
        bool dirichletcond)
  : fluid_(fluiddis, xfluiddis, *solver, *params, *output, isale),
    fluiddis_(fluiddis),
    solver_(solver),
    params_(params),
    output_(output)
{

  interface_.Setup(*fluiddis);
  fluid_.SetSurfaceSplitter(&interface_);

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint
  const Teuchos::RCP<const LINALG::MapExtractor> dbcmaps = fluid_.DirichMaps();
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(interface_.OtherMap());
  maps.push_back(dbcmaps->OtherMap());
  innervelmap_ = LINALG::MultiMapExtractor::IntersectMaps(maps);

  if (dirichletcond)
  {
    // mark all interface velocities as dirichlet values
    fluid_.AddDirichCond(interface_.FSICondMap());
  }

  vector<string> conditions_to_copy;
  conditions_to_copy.push_back("XFEMCoupling");
  fluidxfluidboundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(fluiddis_, "XFEMCoupling", "boundary", "BELE3", conditions_to_copy);
  if (fluidxfluidboundarydis_->NumGlobalNodes() == 0)
  {
    cout << "Empty fluidxfluidboundary discretization detected!" << endl;
  }

  // Note: The nodal structure doesn't change, so we'll create the nodal maps in the constructor.
  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map ffnoderowmap = *fluidxfluidboundarydis_->NodeRowMap();
  const Epetra_Map ffelemrowmap = *fluidxfluidboundarydis_->ElementRowMap();

  // put all boundary nodes and elements onto all processors
  // Create an allreduced Epetra_Map from the given Epetra_Map and give it to all processors
  const Epetra_Map ffnodecolmap = *LINALG::AllreduceEMap(ffnoderowmap);
  const Epetra_Map ffelemcolmap = *LINALG::AllreduceEMap(ffelemrowmap);

  // redistribute nodes and elements to column (ghost) map
  fluidxfluidboundarydis_->ExportColumnNodes(ffnodecolmap);
  fluidxfluidboundarydis_->ExportColumnElements(ffelemcolmap);

  const int err =  fluidxfluidboundarydis_->FillComplete();
  if (err) dserror("FillComplete() returned err=%d",err);

  RCP<Epetra_Map> newcolnodemap = DRT::UTILS::ComputeNodeColMap(fluiddis_, fluidxfluidboundarydis_);
  fluiddis_->Redistribute(*(fluiddis_->NodeRowMap()), *newcolnodemap);

  DRT::UTILS::PrintParallelDistribution(*fluidxfluidboundarydis_);

  fluidxfluidboundaryoutput_ = rcp(new IO::DiscretizationWriter(fluidxfluidboundarydis_));

  // create fluid-xfluid-interface DOF vectors
  fxfidispnp_    = LINALG::CreateVector(*fluidxfluidboundarydis_->DofRowMap(),true);
  fxfivelnp_     = LINALG::CreateVector(*fluidxfluidboundarydis_->DofRowMap(),true); // physical fluid velocity
  fxfitrueresnp_ = LINALG::CreateVector(*fluidxfluidboundarydis_->DofRowMap(),true);
  fxfidispn_   = LINALG::CreateVector(*fluidxfluidboundarydis_->DofRowMap(),true);
  fxfiveln_    = LINALG::CreateVector(*fluidxfluidboundarydis_->DofRowMap(),true);
  fxfivelnm_   = LINALG::CreateVector(*fluidxfluidboundarydis_->DofRowMap(),true);
  fxfiaccnp_   = LINALG::CreateVector(*fluidxfluidboundarydis_->DofRowMap(),true);
  fxfiaccn_    = LINALG::CreateVector(*fluidxfluidboundarydis_->DofRowMap(),true);

  PrepareFluidXFluidBoundaryDis();
  fluid_.PrepareFluidXFluidBoundaryDofset(fluidxfluidboundarydis_);
  fluid_.PrepareFluidXFluidBoundaryDis(fluidxfluidboundarydis_);
  fluid_.PrepareTimeLoop(fluidxfluidboundarydis_);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFluidImpl::InitialGuess()
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFluidImpl::RHS()
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFluidImpl::TrueResidual()
{
  return fluid_.TrueResidual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFluidImpl::Velnp()
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFluidImpl::Velaf()
{
  return fluid_.Velaf();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFluidImpl::Veln()
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFluidImpl::Dispnp()
{
  return fluid_.Dispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFluidImpl::ConvectiveVel()
{
  if (fluid_.GridVel() == Teuchos::null)
    return fluid_.Velnp(); // no moving mesh present
  else
  {
    // make an intermediate copy of velnp
    Teuchos::RCP<Epetra_Vector> convel = rcp(new Epetra_Vector(*(fluid_.Velnp())));
    // now subtract the grid velocity
    convel->Update(-1.0,*(fluid_.GridVel()),1.0);

    return convel;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidXFluidImpl::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MapExtractor> ADAPTER::FluidXFluidImpl::GetDBCMapExtractor()
{
  return fluid_.DirichMaps();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::TimeLoop()
{
  fluid_.Integrate(fluidxfluidboundarydis_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::Update()
{
  fluid_.TimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::StatisticsAndOutput()
{
  fluid_.StatisticsAndOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::Output()
{
  fluid_.Output();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::NonlinearSolve()
{
  cout << "ADAPTER::FluidXFluidImpl::NonlinearSolve()" << endl;
  fluid_.NonlinearSolve(fluidxfluidboundarydis_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidXFluidImpl::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidXFluidImpl::VelocityRowMap()
{
  return fluid_.VelocityRowMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidXFluidImpl::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_.Setup(*fluiddis_->DofRowMap(),mm,LINALG::SplitMap(*fluiddis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidXFluidImpl::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidXFluidImpl::TimeScaling() const
{
  if (params_->get<bool>("interface second order"))
  {
    return 2./fluid_.Dt();
  }
  else
    return 1./fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidXFluidImpl::Time() const
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidXFluidImpl::Step() const
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidXFluidImpl::Dt() const
{
  return fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFluidImpl::ExtractInterfaceForces()
{
  return interface_.ExtractFSICondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFluidImpl::ExtractInterfaceFluidVelocity()
{
  return interface_.ExtractFSICondVector(fluid_.Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFluidImpl::ExtractInterfaceVeln()
{
  return interface_.ExtractFSICondVector(fluid_.Veln());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertFSICondVector(ivel,fluid_.Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSICondVector(Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = TimeScaling();
  fcx->Update(-timescale*fluid_.Dt(),*veln,timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSICondVector(Veln());

  // We convert Delta u(n+1,i+1) to Delta d(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = 1./TimeScaling();
  fcx->Update(fluid_.Dt(),*veln,timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidXFluidImpl::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFluidImpl::IntegrateInterfaceShape()
{
  return interface_.ExtractFSICondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidXFluidImpl::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::FluidXFluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFluidImpl::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
{
  return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}


///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno)
{
  fluid_.SetInitialFlowField(fluidxfluidboundarydis_,initfield,startfuncno);
  return;
}


///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
void ADAPTER::FluidXFluidImpl::PrepareFluidXFluidBoundaryDis()
{
  // put vectors into boundary discretization (SetState generates col vector automatically)
  fluidxfluidboundarydis_->SetState("idispcolnp",fxfidispnp_);
  fluidxfluidboundarydis_->SetState("idispcoln" ,fxfidispn_);
  fluidxfluidboundarydis_->SetState("ivelcolnp" ,fxfivelnp_);
  fluidxfluidboundarydis_->SetState("ivelcoln"  ,fxfiveln_);
  fluidxfluidboundarydis_->SetState("ivelcolnm" ,fxfivelnm_);
  fluidxfluidboundarydis_->SetState("iacccoln"  ,fxfiaccn_);
}


#endif
