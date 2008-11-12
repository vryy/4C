/*----------------------------------------------------------------------*/
/*!
\file passive_scatra_algorithm.cpp

\brief Transport of passive scalars in Navier-Stokes velocity field

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "passive_scatra_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::PassiveScaTraAlgorithm::PassiveScaTraAlgorithm(Epetra_Comm& comm, const Teuchos::ParameterList& prbdyn)
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn)
{
 // no stuff to add here at the moment
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::PassiveScaTraAlgorithm::~PassiveScaTraAlgorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::TimeLoop()
{
  // write out inital state
  // Output();

  // time loop (no-subcycling at the moment)
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    // solve nonlinear Navier-Stokes system
    DoFluidStep();

    // solve transport (convection-diffusion) equations for passive scalar 
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
void SCATRA::PassiveScaTraAlgorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n   FLUID SOLVER  \n******************\n";
  }

  FluidField().PrepareTimeStep();

  // transfer the initial(!!) convective velocity
  //(fluid initial field was set inside the constructor of fluid base class)
  if (Step()==1)
    ScaTraField().SetVelocityField(2,ConvectiveVelocity());

  // prepare time step (+ initialize one-step-theta scheme correctly with 
  // velocity given above)
  ScaTraField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::DoFluidStep()
{
  // solve nonlinear Navier-Stokes system
  FluidField().NonlinearSolve();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::DoTransportStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n TRANSPORT SOLVER \n******************\n";
  }
  
  // get new velocity from Navier-Stokes solver
  GetCurrentFluidVelocity();

  // transfer convective velocity to scalar transport field solver
  ScaTraField().SetVelocityField(2,ConvectiveVelocity());

  // solve the linear convection-diffusion equation(s)
  ScaTraField().Solve();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::Update()
{
  FluidField().Update();
  ScaTraField().Update();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::PassiveScaTraAlgorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField().Output();
  FluidField().LiftDrag();
  ScaTraField().Output();

  return;
}


#endif // CCADISCRET
