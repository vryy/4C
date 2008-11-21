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
    RCP<DRT::Discretization>      fluiddis,
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn),
   fluiddiscret_(fluiddis)
{
  // maximum number of iterations and tolerance for outer iteration
  itmax_ = prbdyn.get<int>   ("ITEMAX");
  ittol_ = prbdyn.get<double>("CONVTOL");

  // flag for printing out mean values of temperature and density
  outmean_ = prbdyn.get<string>("OUTMEAN");

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
  if (Comm().MyPID()==0)
    cout<<"\n**********************\n STATIONARY LOW-MACH-NUMBER FLOW SOLVER \n**********************\n";

  // prepare time step
  PrepareTimeStep();

  // do outer iteration loop
  OuterLoop();

  // write output to screen and files
  Output();

return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::TimeLoop()
{
  // do initial calculations
  InitialCalculations();

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
void LOMA::Algorithm::InitialCalculations()
{
  // compute initial density field using initial temperature + therm. pressure
  // temperature stored at phinp in SCATRA is used -> densnp is set
  ScaTraField().ComputeDensity();

  // initially set density at 0 (i.e., densn)
  // furthermore, set density at -1 (i.e., densnm) for BDF2 (zero vector)
  ScaTraField().UpdateDensity();

  // get initial velocity (and pressure) field
  // use value at n+1 for all time-integration schemes, i.e., also for
  // generalized-alpha time integration
  GetFluidVelPressNp();

  // compute initial convective density-weighted velocity field for scalar
  // transport solver using initial fluid velocity (and pressure) field
  // for generalized-alpha time-integration scheme, velocity at n+1
  // is weighted by density at n+alpha_F, which is identical to density
  // at n+1, since density at n was set equal to n+1 above
  ScaTraField().SetLomaVelocity(VelocityPressureNp(),fluiddiscret_);

  // write initial fields
  Output();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::PrepareTimeStep()
{
  // predict density field (from second time step on)
  if (Step()> 1) ScaTraField().PredictDensity();

  // prepare temperature time step (+ initialize one-step-theta and
  // generalized-alpha schemes correctly by computing initial time
  // derivatives of temperature dphi/dt_(0))
  ScaTraField().PrepareTimeStep();

  // get density at n+1 and n
  GetDensityNp();
  GetDensityN();
  if (ScaTraField().MethodName()==INPAR::SCATRA::timeint_gen_alpha)
  {
    // get density time derivative at n
    GetDensityDtN();

     // set density at n+1 and n as well as density time derivative at n
    FluidField().SetTimeLomaFields(DensityNp(),DensityN(),DensityDtN());
  }
  else
  {
    // get density at n-1
    GetDensityNm();

     // set density at n+1 , n and n-1
    FluidField().SetTimeLomaFields(DensityNp(),DensityN(),DensityNm());
  }

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

    if (ScaTraField().MethodName()==INPAR::SCATRA::timeint_gen_alpha)
    {
      // get velocity (and pressure) field at intermediate time step n+alpha_F
      GetFluidVelPressAf();

      // set field vectors: density-weighted convective velocity + density
      ScaTraField().SetLomaVelocity(VelocityPressureAf(),fluiddiscret_);
    }
    else
    {
      // get velocity (and pressure) field at time step n+1
      GetFluidVelPressNp();

      // set field vectors: density-weighted convective velocity + density
      ScaTraField().SetLomaVelocity(VelocityPressureNp(),fluiddiscret_);
    }

    // solve transport equation for temperature
    if (Comm().MyPID()==0) cout<<"\n**********************\n  TEMPERATURE SOLVER \n**********************\n";
    ScaTraField().Solve();

    // compute density using current temperature + thermodynamic pressure
    ScaTraField().ComputeDensity();

    // get current density at n+1
    GetDensityNp();

    if (ScaTraField().MethodName()==INPAR::SCATRA::timeint_gen_alpha)
    {
      // compute time derivative of density
      ScaTraField().ComputeDensityDerivative();

      // get density time derivative at n+1
      GetDensityDtNp();

       // set density at n+1 and n as well as density time derivative at n
      FluidField().SetGenAlphaIterLomaFields(DensityNp(),DensityDtNp());
    }
    else
    {
      // set current density
      FluidField().SetIterLomaFields(DensityNp());
    }

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0) cout<<"\n*********************\n     FLOW SOLVER \n*********************\n";
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
  // update density
  ScaTraField().UpdateDensity();

  // update temperature
  ScaTraField().Update();

  // get density at n+1 and n
  GetDensityNp();
  GetDensityN();
  if (ScaTraField().MethodName()==INPAR::SCATRA::timeint_gen_alpha)
  {
    // get density time derivative at n
    GetDensityDtN();

     // set density at n+1 and n as well as density time derivative at n
    FluidField().SetTimeLomaFields(DensityNp(),DensityN(),DensityDtN());
  }
  else
  {
    // get density at n-1
    GetDensityNm();

     // set density at n+1 , n and n-1
    FluidField().SetTimeLomaFields(DensityNp(),DensityN(),DensityNm());
  }

  // update fluid
  FluidField().Update();

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
  FluidField().StatisticsAndOutput();

  ScaTraField().Output();
  if (outmean_=="Yes") ScaTraField().OutputMeanTempAndDens();

  return;
}


#endif // CCADISCRET
