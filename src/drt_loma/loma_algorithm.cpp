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
  // check whether stationary algorithm was chosen
  if (FluidField().TimIntScheme() == timeint_stationary)
    dserror("no stationary algorithm supported currently");

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
  // set initial velocity field for scalar transport solver
  ScaTraField().SetVelocityField(FluidField().Velnp(),
                                 FluidField().SgVelVisc(),
                                 FluidField().Discretization());

  // set initial value of thermodynamic pressure in SCATRA
  ScaTraField().SetInitialThermPressure();

  // compute initial mass, initial thermodynamic pressure (only in case of mass
  // conservation) and initial time derivative of thermodynamic pressure
  if (consthermpress_=="No_energy")
    ScaTraField().ComputeInitialThermPressureDeriv();
  else if (consthermpress_=="No_mass")
  {
    ScaTraField().ComputeInitialMass();
    ScaTraField().ComputeThermPressureFromMassCons();
  }

  // write initial fields
  //Output();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LOMA::Algorithm::PrepareTimeStep()
{
  // prepare scalar transport time step (+ initialize one-step-theta and
  // generalized-alpha schemes correctly by computing initial time
  // derivatives of scalar dphi/dt_(0) in first time step)
  ScaTraField().PrepareTimeStep();

  // predict thermodynamic pressure and time derivative
  // (if not constant or based on mass conservation)
  if (consthermpress_=="No_energy") ScaTraField().PredictThermPressure();

  // set scalar and thermodynamic pressure at n+1 and n as well as their
  // time derivatives at n, SCATRA trueresidual and number of scalars
  FluidField().SetTimeLomaFields(ScaTraField().Phinp(),
                                 ScaTraField().Phin(),
                                 ScaTraField().PhiDtn(),
                                 ScaTraField().ThermPressNp(),
                                 ScaTraField().ThermPressN(),
                                 ScaTraField().ThermPressDtN(),
                                 ScaTraField().TrueResidual(),
                                 ScaTraField().NumScal());

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
    cout<<"\n******************************************\n  OUTER ITERATION LOOP\n******************************************\n";

  // maximum number of iterations tolerance for outer iteration
  // currently default for turbulent channel flow: only one iteration before sampling
  if (special_flow_ == "loma_channel_flow_of_height_2" && Step() < samstart_ )
       itmax_ = 1;
  else itmax_ = itmaxpre_;

  // evaluate fluid predictor step (currently not performed)
  //FluidField().Predictor();

  while (stopnonliniter==false)
  {
    itnum++;

    // store scalar (and velocity) of previous iteration for convergence check
    ScaTraField().ScalIncNp()->Update(1.0,*ScaTraField().Phinp(),0.0);
    //ScaTraField().VelIncNp()->Update(1.0,*ConVel(),0.0);

    // compute values at intermediate time steps (only for generalized-alpha
    // time-integration scheme) and set respective field vectors for velocity,
    // subgrid velocity/viscosity and discretization
    if (FluidField().TimIntScheme() == timeint_afgenalpha)
    {
      ScaTraField().ComputeIntermediateValues();

      ScaTraField().SetVelocityField(FluidField().Velaf(),
                                     FluidField().SgVelVisc(),
                                     FluidField().Discretization());
    }
    else
      ScaTraField().SetVelocityField(FluidField().Velnp(),
                                     FluidField().SgVelVisc(),
                                     FluidField().Discretization());

    // solve transport equation for temperature
    if (Comm().MyPID()==0) cout<<"\n******************************************\n   TEMPERATURE SOLVER\n******************************************\n";
    ScaTraField().Solve();

    // in case of non-constant thermodynamic pressure: compute
    if (consthermpress_=="No_energy")
      ScaTraField().ComputeThermPressure();
    else if (consthermpress_=="No_mass")
      ScaTraField().ComputeThermPressureFromMassCons();

    // set scalar and thermodynamic pressure values as well as time derivatives
    // at n+1 and number of scalars
    FluidField().SetIterLomaFields(ScaTraField().Phinp(),
                                   ScaTraField().PhiDtnp(),
                                   ScaTraField().ThermPressNp(),
                                   ScaTraField().ThermPressDtNp(),
                                   ScaTraField().NumScal());

    // solve low-Mach-number flow equations
    if (Comm().MyPID()==0) cout<<"\n******************************************\n      FLOW SOLVER\n******************************************\n";
    FluidField().MultiCorrector();

    // check convergence of temperature field
    stopnonliniter = ScaTraField().LomaConvergenceCheck(itnum,itmax_,ittol_);
  }

  // compute values at intermediate time steps (only for generalized-alpha
  // time-integration scheme) and set respective field vectors for velocity,
  // subgrid velocity/viscosity and discretization
  if (FluidField().TimIntScheme() == timeint_afgenalpha)
  {
    ScaTraField().ComputeIntermediateValues();

    ScaTraField().SetVelocityField(FluidField().Velaf(),
                                   FluidField().SgVelVisc(),
                                   FluidField().Discretization());
  }
  else
    ScaTraField().SetVelocityField(FluidField().Velnp(),
                                   FluidField().SgVelVisc(),
                                   FluidField().Discretization());

  // solve transport equation for temperature
  if (Comm().MyPID()==0) cout<<"\n******************************************\n   TEMPERATURE SOLVER\n******************************************\n";
  ScaTraField().Solve();

  // in case of non-constant thermodynamic pressure: compute
  if (consthermpress_=="No_energy")
    ScaTraField().ComputeThermPressure();
  else if (consthermpress_=="No_mass")
    ScaTraField().ComputeThermPressureFromMassCons();

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
  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField().StatisticsAndOutput();

  ScaTraField().Output();

  return;
}


#endif // CCADISCRET
