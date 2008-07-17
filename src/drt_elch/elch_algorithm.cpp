/*----------------------------------------------------------------------*/
/*!
\file elch_algorithm.cpp

\brief Basis of all ELCH algorithms

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "elch_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::Algorithm(Epetra_Comm& comm, ParameterList& prbdyn)
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn)
{
 // no ELCH specific stuff at the moment
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

    // solve transport equations for ion concentrations and electric potential
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
  GetCurrentFluidVelocity();

  // prepare time step
  ScaTraField().PrepareTimeStep();
  // transfer convective velocity
  ScaTraField().SetVelocityField(2,ConvectiveVelocity());
  //ConDifField().SetVelocityField(1,1);

  // solve coupled transport equations for ion concentrations and electric 
  // potential
  ScaTraField().Solve();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::Update()
{
  FluidField().Update();
  ScaTraField().Update();
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
  ScaTraField().Output();

  // debug IO
#if 0
  // print out convective velocity vector
  cout<<*velocitynp_<<endl;
#endif

  return;
}


#endif // CCADISCRET
