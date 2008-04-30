
#ifdef CCADISCRET

#include "elch_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::Algorithm(Epetra_Comm& comm)
:  FluidBaseAlgorithm(DRT::Problem::Instance()->FluidDynamicParams(),false),
   ConDifBaseAlgorithm(DRT::Problem::Instance()->FluidDynamicParams()),
   comm_(comm),
   step_(0),
   time_(0.0)
{
  // taking time loop control parameters out of fluid dynamics section
  const Teuchos::ParameterList& fluiddyn = DRT::Problem::Instance()->FluidDynamicParams();
  // maximum simulation time
  maxtime_=fluiddyn.get<double>("MAXTIME");
  // maximum number of timesteps
  nstep_ = fluiddyn.get<int>("NUMSTEP");
  // time step size
  dt_ = fluiddyn.get<double>("TIMESTEP");

  // get RCP to actual velocity field (time n+1)
  velocitynp_=FluidField().ExtractVelocityPart(FluidField().Velnp());

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::~Algorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::TimeLoop()
{
    // solve con-dif equation without coupling to Navier-Stokes
    // ConDifField().SetVelocityField(1,1);
    // ConDifField().Integrate();

  // time loop
  while (NotFinished())
  {
    // prepare next time step (only ELCH and fluid)
    PrepareTimeStep();

    // solve nonlinear Navier-Stokes system
    DoFluidStep();

    // solve convection-diffusion equation
    DoTransportStep();

    // update all field solvers
    Update();

    // write output to screen and files
    Output();

  } // time loop

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n   FLUID SOLVER  \n******************\n";
  }

  FluidField().PrepareTimeStep();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::DoFluidStep()
{
  // solve nonlinear Navier-Stokes system
  FluidField().NonlinearSolve();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::DoTransportStep()
{
  
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n TRANSPORT SOLVER \n******************\n";
  }
  
  // get new velocity from Navier-Stokes solver
  velocitynp_=FluidField().ExtractVelocityPart(FluidField().Velnp());

  // prepare time step
  ConDifField().PrepareTimeStep();
  // transfer convective velocity
  //ConDifField().SetVelocityField(2,velocitynp_);
  ConDifField().SetVelocityField(1,1);

  // solve convection-diffusion equation
  ConDifField().Solve();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::Update()
{
  FluidField().Update();
  ConDifField().Update();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField().Output();
  FluidField().LiftDrag();
  ConDifField().Output();

  // debug IO
#if 0
  // print out convective velocity vector
  cout<<*velocitynp_<<endl;
#endif

  return;
}


#endif // CCADISCRET
