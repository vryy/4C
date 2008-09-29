/*----------------------------------------------------------------------*/
/*!
\file loma_algorithm.cpp

\brief Basis of all LOMA algorithms

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "loma_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LOMA::Algorithm::Algorithm(
    Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn)
{
  // maximum number of iterations and tolerance for outer iteration
  itmax_ = prbdyn.get<int>   ("ITEMAX");
  ittol_ = prbdyn.get<double>("CONVTOL");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LOMA::Algorithm::~Algorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::SolveStationaryProblem()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
    cout<<"\n**********************\n STATIONARY OUTER ITERATION LOOP \n**********************\n";

  // set boundary conditions for scalar transport solver
  ScaTraField().PrepareTimeStep();

  // set boundary conditions for fluid solver
  FluidField().PrepareTimeStep();

  while (stopnonliniter==false)
  {
    itnum++;

    // get new velocity field
    GetCurrentFluidVelocity();

    // set field vectors: density-weighted convective velocity + density
    ScaTraField().SetLomaVelocity(ConvectiveVelocity());

    // solve transport equation for temperature
    if (Comm().MyPID()==0)
      cout<<"\n**********************\n  TEMPERATURE SOLVER \n**********************\n";
    ScaTraField().Solve();

    // compute density using current temperature + thermodynamic pressure
    ScaTraField().ComputeDensity();

    // get temperature
    GetCurrentDensity();

    // set field vectors: density
    FluidField().SetIterLomaFields(Density());

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0)
      cout<<"\n*********************\n     FLOW SOLVER \n*********************\n";
    FluidField().NonlinearSolve();

    // check convergence of density
    stopnonliniter = ScaTraField().DensityConvergenceCheck(itnum,itmax_,ittol_);
  }

  Output();

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::TimeLoop()
{
  // compute initial density field using initial temperature + therm. pressure
  ScaTraField().ComputeDensity();

  // initially set density at n and n-1
  ScaTraField().UpdateDensity();

  // get initial velocity field
  GetCurrentFluidVelocity();

  // compute initial convective density-weighted velocity field for scalar 
  // transport solver using initial fluid velocity field
  ScaTraField().SetLomaVelocity(ConvectiveVelocity());

  // time loop
  while (NotFinished())
  {
    IncrementTimeAndStep();

    // prepare time step
    PrepareTimeStep();

    // do outer iteration loop
    OuterLoop();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

  } // time loop

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::PrepareTimeStep()
{
  // predict density field (from second time step on)
  if (Step()> 1) ScaTraField().PredictDensity();

  // prepare temperature time step (+ initialize one-step-theta scheme correctly)
  ScaTraField().PrepareTimeStep();

  // get temperature at n+1, n and n-1
  GetCurrentDensity();
  GetNDensity();
  GetNmDensity();

  // set density at n+1, n and n-1
  FluidField().SetTimeLomaFields(Density(),NDensity(),NmDensity());

  // prepare fluid time step, particularly predict velocity field
  FluidField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
    cout<<"\n**********************\n OUTER ITERATION LOOP \n**********************\n";

  while (stopnonliniter==false)
  {
    itnum++;

    // get new velocity field
    GetCurrentFluidVelocity();

    // set field vectors: density-weighted convective velocity + density
    ScaTraField().SetLomaVelocity(ConvectiveVelocity());

    // solve transport equation for temperature
    if (Comm().MyPID()==0)
      cout<<"\n**********************\n  TEMPERATURE SOLVER \n**********************\n";
    ScaTraField().Solve();

    // compute density using current temperature + thermodynamic pressure
    ScaTraField().ComputeDensity();

    // get temperature
    GetCurrentDensity();

    // set field vectors: density
    FluidField().SetIterLomaFields(Density());

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0)
      cout<<"\n*********************\n     FLOW SOLVER \n*********************\n";
    FluidField().NonlinearSolve();

    // check convergence of density
    stopnonliniter = ScaTraField().DensityConvergenceCheck(itnum,itmax_,ittol_);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::Update()
{
  FluidField().Update();
  ScaTraField().Update();
  ScaTraField().UpdateDensity();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::Output()
{
  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField().Output();
  FluidField().LiftDrag();
  ScaTraField().Output();

  return;
}


#endif // CCADISCRET
