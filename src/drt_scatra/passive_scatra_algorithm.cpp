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
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PassiveScaTraAlgorithm::PassiveScaTraAlgorithm(Epetra_Comm& comm, const Teuchos::ParameterList& prbdyn)
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn)
{
 // no stuff to add here at the moment
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PassiveScaTraAlgorithm::~PassiveScaTraAlgorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PassiveScaTraAlgorithm::PassiveScaTraAlgorithm::TimeLoop()
{
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
void PassiveScaTraAlgorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n   FLUID SOLVER  \n******************\n";
  }

  FluidField().PrepareTimeStep();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PassiveScaTraAlgorithm::DoFluidStep()
{
  // solve nonlinear Navier-Stokes system
  FluidField().NonlinearSolve();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PassiveScaTraAlgorithm::DoTransportStep()
{
  
  if (Comm().MyPID()==0)
  {
    cout<<"\n******************\n TRANSPORT SOLVER \n******************\n";
  }
  
  // get new velocity from Navier-Stokes solver
  GetCurrentFluidVelocity();

  // prepare time step
  ScaTraField().PrepareTimeStep();
  // transfer convective velocity to scalar transport field solver
  ScaTraField().SetVelocityField(2,ConvectiveVelocity());

  // solve the linear convection-diffusion equation(s)
  ScaTraField().Solve();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PassiveScaTraAlgorithm::Update()
{
  FluidField().Update();
  ScaTraField().Update();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PassiveScaTraAlgorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField().Output();
  FluidField().LiftDrag();
  ScaTraField().Output();

  // debug IO
#if 0
  // print out convective velocity vector
  cout<<*velocitynp_<<endl;
#endif

  return;
}


#endif // CCADISCRET
