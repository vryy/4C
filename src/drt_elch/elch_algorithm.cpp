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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::Algorithm(
    Epetra_Comm& comm, 
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn),
   outmean_(Teuchos::getIntegralValue<int>(prbdyn,"OUTMEAN"))
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::~Algorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::TimeLoop()
{
  // write out inital state
  // Output();

  // compute error for problems with analytical solution
  ScaTraField().EvaluateErrorComparedToAnalyticalSol();

  // time loop
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    // solve nonlinear Navier-Stokes system
    DoFluidStep();

    // solve transport equations for ion concentrations and electric potential
    DoTransportStep();

    // update all single field solvers
    Update();

    // compute error for problems with analytical solution
    ScaTraField().EvaluateErrorComparedToAnalyticalSol();

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

  // transfer convective velocity
  ScaTraField().SetVelocityField(2,ConvectiveVelocity());
  //ScaTraField().SetVelocityField(1,1);

  // solve coupled transport equations for ion concentrations and electric 
  // potential
  ScaTraField().NonlinearSolve();
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
  if (outmean_) 
    ScaTraField().OutputMeanTempAndDens();

  return;
}


#endif // CCADISCRET
