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
    Epetra_Comm&                  comm,
    const Teuchos::ParameterList& prbdyn
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false)
{
  // (preliminary) maximum number of iterations and tolerance for outer iteration
  ittol_    = prbdyn.get<double>("CONVTOL");
  itmaxpre_ = prbdyn.get<int>("ITEMAX");

  // flag for constant thermodynamic pressure
  consthermpress_ = prbdyn.get<string>("CONSTHERMPRESS");

  // flag for special flow and start of sampling period from fluid parameter list
  special_flow_ = prbdyn.get<string>("CANONICAL_FLOW");
  samstart_     = prbdyn.get<int>("SAMPLING_START");

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
/*void LOMA::Algorithm::SolveStationaryProblem()
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
}*/


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

    // do outer iteration loop for particular type of algorithm
    if (FluidField().TimIntScheme() == timeint_afgenalpha)
      GenAlphaOuterLoop();
    else if (FluidField().TimIntScheme() == timeint_one_step_theta or
             FluidField().TimIntScheme() == timeint_bdf2)
      OSTBDF2OuterLoop();
    else dserror("desired type of low-Mach-number algorithm not supported");

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
  // set initial velocity field for evaluation of initial scalar time
  // derivative in SCATRA
  ScaTraField().SetVelocityField(FluidField().Velnp(),
                                 Teuchos::null,
                                 Teuchos::null,
                                 FluidField().Discretization());

  // set initial value of thermodynamic pressure in SCATRA
  ScaTraField().SetInitialThermPressure();

  // energy conservation: compute initial time derivative of therm. pressure
  // mass conservation: compute initial mass (initial time deriv. assumed zero)
  // constant thermodynamic pressure: for generalized-alpha
  // time-integration scheme, values at n+alpha_F and n+alpha_M are set
  // once and for all simulation time
  if (consthermpress_=="No_energy")
    ScaTraField().ComputeInitialThermPressureDeriv();
  else if (consthermpress_=="No_mass")
    ScaTraField().ComputeInitialMass();
  else
    ScaTraField().ComputeThermPressureIntermediateValues();

  // set initial scalar field and thermodynamic pressure for evaluation of
  // Neumann boundary conditions in FLUID at beginning of first time step
  FluidField().SetTimeLomaFields(ScaTraField().Phinp(),
                                 ScaTraField().ThermPressNp(),
                                 null,
                                 ScaTraField().NumScal());

  // write initial fields
  //Output();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::PrepareTimeStep()
{
  // prepare scalar transport time step
  // (+ computation of initial scalar time derivative in first time step)
  ScaTraField().PrepareTimeStep();

  // predict thermodynamic pressure and time derivative
  // (if not constant or based on mass conservation)
  if (consthermpress_=="No_energy") ScaTraField().PredictThermPressure();

  // prepare fluid time step, among other things, predict velocity field
  FluidField().PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::GenAlphaOuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0) cout<<"\n************************\n  OUTER ITERATION LOOP\n************************\n";

  // maximum number of iterations tolerance for outer iteration
  // currently default for turbulent channel flow: only one iteration before sampling
  if (special_flow_ == "loma_channel_flow_of_height_2" && Step() < samstart_ )
       itmax_ = 1;
  else itmax_ = itmaxpre_;

  // evaluate fluid predictor step (currently not performed)
  //FluidField().Predictor();

  // set respective field vectors for velocity/pressure, acceleration
  // and discretization
  ScaTraField().SetVelocityField(FluidField().Velaf(),
                                 FluidField().Accam(),
                                 Teuchos::null,
                                 FluidField().Discretization());

  // compute scalar values at intermediate time steps
  ScaTraField().ComputeIntermediateValues();

  // initially solve scalar transport equation
  if (Comm().MyPID()==0) cout<<"\n************************\n     SCALAR SOLVER\n************************\n";
  ScaTraField().Solve();

  // update scalar time derivative after first SCATRA solution
  ScaTraField().UpdateTimeDerivative();

  while (stopnonliniter==false)
  {
    itnum++;

    // store scalar from first solution for convergence check
    ScaTraField().ScalIncNp()->Update(1.0,*ScaTraField().Phinp(),0.0);

    // compute scalar values at intermediate time steps
    ScaTraField().ComputeIntermediateValues();

    // in case of non-constant thermodynamic pressure: compute and update
    if (consthermpress_=="No_energy")
    {
      ScaTraField().ComputeThermPressure();

      // update time derivative of thermodynamic pressure after solution
      ScaTraField().UpdateThermPressureTimeDerivative();

      // compute values of therm. pressure at intermediate time steps
      ScaTraField().ComputeThermPressureIntermediateValues();
    }
    else if (consthermpress_=="No_mass")
    {
      ScaTraField().ComputeThermPressureFromMassCons();

      // update time derivative of thermodynamic pressure after solution
      ScaTraField().UpdateThermPressureTimeDerivative();

      // compute values of therm. pressure at intermediate time steps
      ScaTraField().ComputeThermPressureIntermediateValues();
    }

    // set scalar and thermodynamic pressure values as well as time derivatives
    // at n+alpha_F and n+alpha_M, respectively, and number of scalars
    FluidField().SetIterLomaFields(ScaTraField().Phiaf(),
                                   ScaTraField().Phiam(),
                                   ScaTraField().Phidtam(),
                                   ScaTraField().ThermPressAf(),
                                   ScaTraField().ThermPressAm(),
                                   ScaTraField().ThermPressDtAm(),
                                   ScaTraField().NumScal());

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0) cout<<"\n************************\n      FLOW SOLVER\n************************\n";
    FluidField().MultiCorrector();

    // set respective field vectors for velocity/pressure, acceleration
    // and discretization
    ScaTraField().SetVelocityField(FluidField().Velaf(),
                                   FluidField().Accam(),
                                   Teuchos::null,
                                   FluidField().Discretization());

    // solve scalar transport equation
    if (Comm().MyPID()==0) cout<<"\n************************\n     SCALAR SOLVER\n************************\n";
    ScaTraField().Solve();

    // update scalar time derivative after SCATRA solution
    ScaTraField().UpdateTimeDerivative();

    // check convergence of temperature field
    stopnonliniter = ScaTraField().LomaConvergenceCheck(itnum,itmax_,ittol_);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::OSTBDF2OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0) cout<<"\n************************\n  OUTER ITERATION LOOP\n************************\n";

  // maximum number of iterations tolerance for outer iteration
  // currently default for turbulent channel flow: only one iteration before sampling
  if (special_flow_ == "loma_channel_flow_of_height_2" && Step() < samstart_ )
       itmax_ = 1;
  else itmax_ = itmaxpre_;

  // evaluate fluid predictor step (currently not performed)
  //FluidField().Predictor();

  // set respective field vectors for velocity/pressure, history vector
  // and discretization
  ScaTraField().SetVelocityField(FluidField().Velnp(),
                                 FluidField().Hist(),
                                 Teuchos::null,
                                 FluidField().Discretization());

  // initially solve scalar transport equation
  if (Comm().MyPID()==0) cout<<"\n************************\n     SCALAR SOLVER\n************************\n";
  ScaTraField().Solve();

  // update scalar time derivative after first SCATRA solution
  ScaTraField().UpdateTimeDerivative();

  while (stopnonliniter==false)
  {
    itnum++;

    // store scalar from first solution for convergence check
    ScaTraField().ScalIncNp()->Update(1.0,*ScaTraField().Phinp(),0.0);

    // in case of non-constant thermodynamic pressure: compute and update
    if (consthermpress_=="No_energy")
    {
      ScaTraField().ComputeThermPressure();

      // update time derivative of thermodynamic pressure after solution
      ScaTraField().UpdateThermPressureTimeDerivative();
    }
    else if (consthermpress_=="No_mass")
    {
      ScaTraField().ComputeThermPressureFromMassCons();

      // update time derivative of thermodynamic pressure after solution
      ScaTraField().UpdateThermPressureTimeDerivative();
    }

    // set scalar and thermodynamic pressure values as well as time derivatives
    // at n+1 and n, respectively, and number of scalars
    FluidField().SetIterLomaFields(ScaTraField().Phinp(),
                                   ScaTraField().Phin(),
                                   ScaTraField().Phidtnp(),
                                   ScaTraField().ThermPressNp(),
                                   ScaTraField().ThermPressN(),
                                   ScaTraField().ThermPressDtNp(),
                                   ScaTraField().NumScal());

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0) cout<<"\n************************\n      FLOW SOLVER\n************************\n";
    FluidField().MultiCorrector();

    // set respective field vectors for velocity/pressure, history vector
    // and discretization
    ScaTraField().SetVelocityField(FluidField().Velnp(),
                                   FluidField().Hist(),
                                   Teuchos::null,
                                   FluidField().Discretization());

    // solve scalar transport equation
    if (Comm().MyPID()==0) cout<<"\n************************\n     SCALAR SOLVER\n************************\n";
    ScaTraField().Solve();

    // update scalar time derivative after SCATRA solution
    ScaTraField().UpdateTimeDerivative();

    // check convergence of temperature field
    stopnonliniter = ScaTraField().LomaConvergenceCheck(itnum,itmax_,ittol_);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::Update()
{
  // update scalar
  ScaTraField().Update();

  // in case of non-constant thermodynamic pressure: update
  if (consthermpress_=="No_energy" or consthermpress_=="No_mass")
    ScaTraField().UpdateThermPressure();

  // update fluid
  FluidField().Update();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::Output()
{
  // set scalar and thermodynamic pressure at n+1 and SCATRA trueresidual
  // for statistical evaluation and evaluation of Neumann boundary
  // conditions at the beginning of the subsequent time step
  FluidField().SetTimeLomaFields(ScaTraField().Phinp(),
                                 ScaTraField().ThermPressNp(),
                                 ScaTraField().TrueResidual(),
                                 ScaTraField().NumScal());

  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField().StatisticsAndOutput();

  ScaTraField().Output();

  return;
}


#endif // CCADISCRET
