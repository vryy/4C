
/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_fluid_impl.cpp


<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_fluid_fluid_impl.H"
#include "../drt_lib/drt_condition_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidFluidImpl::FluidFluidImpl(
        Teuchos::RCP<DRT::Discretization> embfluiddis,
        Teuchos::RCP<DRT::Discretization> bgfluiddis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<ParameterList> params,
        bool isale,
        bool dirichletcond)
  : fluid_( bgfluiddis,embfluiddis, *solver, *params,  isale),
    embfluiddis_(embfluiddis),
    bgfluiddis_(bgfluiddis),
    solver_(solver),
    params_(params)
    //output_(output)
{
  interface_.Setup(*embfluiddis);
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidImpl::InitialGuess()
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidImpl::RHS()
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidImpl::TrueResidual()
{
  return fluid_.TrueResidual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidImpl::Velnp()
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidImpl::Velaf()
{
  return fluid_.Velaf();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidImpl::Veln()
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidImpl::Dispnp()
{
  return fluid_.Dispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidImpl::ConvectiveVel()
{
//   if (fluid_.GridVel() == Teuchos::null)
//     return fluid_.Velnp(); // no moving mesh present
//   else
//   {
//     // make an intermediate copy of velnp
//     Teuchos::RCP<Epetra_Vector> convel = rcp(new Epetra_Vector(*(fluid_.Velnp())));
//     // now subtract the grid velocity
//     convel->Update(-1.0,*(fluid_.GridVel()),1.0);

//     return convel;
//   }
  dserror(" ConvectiveVel() not implemented! ");
  return null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidFluidImpl::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MapExtractor> ADAPTER::FluidFluidImpl::GetDBCMapExtractor()
{
  return fluid_.DirichMaps();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::TimeLoop()
{
  fluid_.IntegrateFluidFluid();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::Update()
{
  fluid_.TimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::StatisticsAndOutput()
{
  fluid_.StatisticsAndOutput();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::Output()
{
  fluid_.Output();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::NonlinearSolve()
{
  fluid_.NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFluidImpl::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFluidImpl::VelocityRowMap()
{
  //return fluid_.VelocityRowMap();
  dserror(" VelocityRowMapx not implemented! ");
  return null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFluidImpl::PressureRowMap()
{
  //return fluid_.PressureRowMap();
  dserror(" pressurerowmap not implemented! ");
  return null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_.Setup(*embfluiddis_->DofRowMap(),mm,LINALG::SplitMap(*embfluiddis_->DofRowMap(),*mm));
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidFluidImpl::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidFluidImpl::TimeScaling() const
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
void ADAPTER::FluidFluidImpl::ReadRestart(int step)
{
//  fluid_.ReadRestart(step);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidFluidImpl::Time() const
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidFluidImpl::Step() const
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidFluidImpl::Dt() const
{
  return fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFluidImpl::ExtractInterfaceForces()
{
  return interface_.ExtractFSICondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFluidImpl::ExtractInterfaceFluidVelocity()
{
  return interface_.ExtractFSICondVector(fluid_.Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFluidImpl::ExtractInterfaceVeln()
{
  return interface_.ExtractFSICondVector(fluid_.Veln());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  cout << "ApplyInterfaceVelocities " << endl;
  interface_.InsertFSICondVector(ivel,fluid_.Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  cout << "ApplyMeshDisplacement " << endl;
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  cout << " DisplacementToVelocity " << endl;
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSICondVector(Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt

  double timescale = TimeScaling();
   fcx->Update(-timescale*fluid_.Dt(),*veln,timescale);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  cout << "VelocityToDisplacement " << endl;
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSICondVector(Veln());

  // We convert Delta u(n+1,i+1) to Delta d(n+1,i+1) here.
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt

  double timescale = 1./TimeScaling();
  fcx->Update(fluid_.Dt(),*veln,timescale);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FluidFluidImpl::Itemax() const
{
  dserror("itermax nicht implementet!");
//  return fluid_.Itemax();
  return null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::SetItemax(int itemax)
{
  dserror("setitermax nicht implementet!");
//  fluid_.SetItemax(itemax);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFluidImpl::IntegrateInterfaceShape()
{
  //return interface_.ExtractFSICondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
  dserror("IntegrateInterfaceShape not implemented!");
  return null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidFluidImpl::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidFluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidFluidImpl::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
{
  cout << "ExtractVelocityPart nicht implemented" << endl;
//  return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
  return null;
}


///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidImpl::SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno)
{
  fluid_.SetInitialFlowField(initfield,startfuncno);
}




#endif
