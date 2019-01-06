/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_ost.cpp
\brief One-Step-Theta time-integration scheme

\level 1

<pre>
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/
#include "scatra_timint_ost.H"

#include "scatra_timint_meshtying_strategy_base.H"
#include "turbulence_hit_scalar_forcing.H"

#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"

#include "../drt_io/io.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::TimIntOneStepTheta(
    Teuchos::RCP<DRT::Discretization> actdis,          //!< discretization
    Teuchos::RCP<LINALG::Solver> solver,               //!< linear solver
    Teuchos::RCP<Teuchos::ParameterList> params,       //!< parameter list
    Teuchos::RCP<Teuchos::ParameterList> extraparams,  //!< supplementary parameter list
    Teuchos::RCP<IO::DiscretizationWriter> output,     //!< output writer
    const int probnum                                  //!< global problem number
    )
    : ScaTraTimIntImpl(actdis, solver, params, extraparams, output, probnum),
      theta_(params_->get<double>("THETA")),
      fsphinp_(Teuchos::null)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting and
  // other kinds of parallel redistribution.
  // This has to be done before all state vectors are initialized.
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Init()
{
  // initialize base class
  ScaTraTimIntImpl::Init();

  return;
}


/*----------------------------------------------------------------------*
 |  setup time integration                                  rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Setup()
{
  // setup base class
  ScaTraTimIntImpl::Setup();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // fine-scale vector at time n+1
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no or
      turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    fsphinp_ = LINALG::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // set element parameters
  // -------------------------------------------------------------------
  // note: - this has to be done before element routines are called
  //       - order is important here: for safety checks in SetElementGeneralParameters(),
  //         we have to know the time-integration parameters
  SetElementTimeParameter();
  SetElementGeneralParameters();
  SetElementTurbulenceParameters();

  // setup krylov
  PrepareKrylovProjection();

  // -------------------------------------------------------------------
  // initialize forcing for homogeneous isotropic turbulence
  // -------------------------------------------------------------------
  // note: this constructor has to be called after the forcing_ vector has
  //       been initialized; this is done in ScaTraTimIntImpl::Init() called before
  if (special_flow_ == "scatra_forced_homogeneous_isotropic_turbulence")
  {
    if (extraparams_->sublist("TURBULENCE MODEL").get<std::string>("SCALAR_FORCING") == "isotropic")
    {
      homisoturb_forcing_ = Teuchos::rcp(new SCATRA::HomIsoTurbScalarForcing(this));
      // initialize forcing algorithm
      homisoturb_forcing_->SetInitialSpectrum(
          DRT::INPUT::IntegralValue<INPAR::SCATRA::InitialField>(*params_, "INITIALFIELD"));
    }
  }

  // -------------------------------------------------------------------
  // initialize multi-scale material if necessary
  // -------------------------------------------------------------------
  if (macro_scale_)
  {
    // create parameter list for macro-scale elements
    Teuchos::ParameterList eleparams;

    // set action
    eleparams.set<int>("action", SCATRA::micro_scale_initialize);

    // loop over macro-scale elements
    discret_->Evaluate(
        eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  }

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    gjb 08/08 |
*-----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::~TimIntOneStepTheta() { return; }


/*----------------------------------------------------------------------*
 |  set time parameter for element evaluation (usual call)   ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetElementTimeParameter(bool forcedincrementalsolver) const
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", SCATRA::set_time_parameter);
  eleparams.set<bool>("using generalized-alpha time integration", false);
  eleparams.set<bool>("using stationary formulation", false);
  if (forcedincrementalsolver == false)
    eleparams.set<bool>("incremental solver", incremental_);
  else
    eleparams.set<bool>("incremental solver", true);

  eleparams.set<double>("time-step length", dta_);
  eleparams.set<double>("total time", time_);
  eleparams.set<double>("time factor", theta_ * dta_);
  eleparams.set<double>("alpha_F", 1.0);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  return;
}


/*--------------------------------------------------------------------------*
 | set time for evaluation of POINT -Neumann boundary conditions   vg 12/08 |
 *--------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetTimeForNeumannEvaluation(Teuchos::ParameterList& params)
{
  params.set("total time", time_);
  return;
}


/*----------------------------------------------------------------------*
| Print information about current time step to screen                   |
*-----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::PrintTimeStepInfo()
{
  if (myrank_ == 0 and not micro_scale_)
  {
    std::cout << std::endl
              << "TIME: " << std::setw(11) << std::setprecision(4) << std::scientific << time_
              << "/" << std::setw(11) << std::setprecision(4) << std::scientific << maxtime_
              << "  DT = " << std::setw(11) << std::setprecision(4) << std::scientific << dta_
              << "  " << MethodTitle() << " (theta = " << std::setw(3) << std::setprecision(2)
              << theta_ << ") STEP = " << std::setw(4) << step_ << "/" << std::setw(4) << stepmax_
              << std::endl;
  }
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetOldPartOfRighthandside()
{
  // call base class routine
  ScaTraTimIntImpl::SetOldPartOfRighthandside();

  // hist_ = phin_ + dt*(1-Theta)*phidtn_
  hist_->Update(1.0, *phin_, dta_ * (1.0 - theta_), *phidtn_, 0.0);

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ExplicitPredictor() const
{
  // call base class routine
  ScaTraTimIntImpl::ExplicitPredictor();

  // predict discrete solution variables
  phinp_->Update(dta_, *phidtn_, 1.0);

  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                             |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AddNeumannToResidual()
{
  residual_->Update(theta_ * dta_, *neumann_loads_, 1.0);
  return;
}


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AVM3Separation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // AVM3 separation
  Sep_->Multiply(false, *phinp_, *fsphinp_);

  // set fine-scale vector
  discret_->SetState("fsphinp", fsphinp_);

  return;
}


/*----------------------------------------------------------------------*
 | dynamic Smagorinsky model                           rasthofer  08/12 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::DynamicComputationOfCs()
{
  if (turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
  {
    // perform filtering and computation of Prt
    // compute averaged values for LkMk and MkMk
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    DynSmag_->ApplyFilterForDynamicComputationOfPrt(
        phinp_, 0.0, dirichtoggle, *extraparams_, nds_vel_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | dynamic Smagorinsky model                           krank  09/13     |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::DynamicComputationOfCv()
{
  if (turbmodel_ == INPAR::FLUID::dynamic_vreman)
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    Vrem_->ApplyFilterForDynamicComputationOfDt(phinp_, 0.0, dirichtoggle, *extraparams_, nds_vel_);
  }

  return;
}


/*--------------------------------------------------------------------------*
 | add global state vectors specific for time-integration scheme   vg 11/08 |
 *--------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver)
{
  // call base class routine
  ScaTraTimIntImpl::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);

  discret_->SetState("hist", hist_);
  discret_->SetState("phinp", phinp_);

  return;
}


/*----------------------------------------------------------------------*
 | compute time derivative                                     vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ComputeTimeDerivative()
{
  // call base class routine
  ScaTraTimIntImpl::ComputeTimeDerivative();

  // time derivative of phi:
  // phidt(n+1) = (phi(n+1)-phi(n)) / (theta*dt) + (1-(1/theta))*phidt(n)
  const double fact = 1.0 / (theta_ * dta_);
  phidtnp_->Update(fact, *phinp_, -fact, *hist_, 0.0);

  // We know the first time derivative on Dirichlet boundaries
  // so we do not need an approximation of these values!
  // However, we do not want to break the linear relationship
  // as stated above. We do not want to set Dirichlet values for
  // dependent values like phidtnp_. This turned out to be inconsistent.
  // ApplyDirichletBC(time_,Teuchos::null,phidtnp_);

  return;
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Update(const int num)
{
  // compute time derivative at time n+1
  ComputeTimeDerivative();

  // call base class routine
  ScaTraTimIntImpl::Update(num);

  // compute flux vector field for later output BEFORE time shift of results
  // is performed below !!
  if (calcflux_domain_ != INPAR::SCATRA::flux_none or
      calcflux_boundary_ != INPAR::SCATRA::flux_none)
  {
    if (DoOutput() or DoBoundaryFluxStatistics()) CalcFlux(true, num);
  }

  // after the next command (time shift of solutions) do NOT call
  // ComputeTimeDerivative() anymore within the current time step!!!

  // solution of this step becomes most recent solution of the last step
  phin_->Update(1.0, *phinp_, 0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  phidtn_->Update(1.0, *phidtnp_, 0.0);

  // call time update of forcing routine
  if (homisoturb_forcing_ != Teuchos::null) homisoturb_forcing_->TimeUpdateForcing();

  // update micro scale in multi-scale simulations if necessary
  if (macro_scale_)
  {
    // create parameter list for macro-scale elements
    Teuchos::ParameterList eleparams;

    // set action
    eleparams.set<int>("action", SCATRA::micro_scale_update);

    // loop over macro-scale elements
    discret_->Evaluate(eleparams);
  }

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::OutputRestart() const
{
  // call base class routine
  ScaTraTimIntImpl::OutputRestart();

  // additional state vectors that are needed for One-Step-Theta restart
  output_->WriteVector("phidtn", phidtn_);
  output_->WriteVector("phin", phin_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input)
{
  // call base class routine
  ScaTraTimIntImpl::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  time_ = reader->ReadDouble("time");
  step_ = reader->ReadInt("step");

  if (myrank_ == 0)
    std::cout << "Reading ScaTra restart data (time=" << time_ << " ; step=" << step_ << ")"
              << std::endl;

  // read state vectors that are needed for One-Step-Theta restart
  reader->ReadVector(phinp_, "phinp");
  reader->ReadVector(phin_, "phin");
  reader->ReadVector(phidtn_, "phidtn");

  ReadRestartProblemSpecific(step, *reader);

  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no or
      turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    AVM3Preparation();

  // read restart on micro scale in multi-scale simulations if necessary
  if (macro_scale_)
  {
    // create parameter list for macro-scale elements
    Teuchos::ParameterList eleparams;

    // set action
    eleparams.set<int>("action", SCATRA::micro_scale_read_restart);

    // loop over macro-scale elements
    discret_->Evaluate(eleparams);
  }

  return;
}


/*-----------------------------------------------------------------------------------------------------------*
 | calculate consistent initial scalar time derivatives in compliance with initial scalar field fang
 09/15 |
 *-----------------------------------------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::CalcInitialTimeDerivative()
{
  // standard general element parameter without stabilization
  SetElementGeneralParameters(true);

  // we also have to modify the time-parameter list (incremental solve)
  // actually we do not need a time integration scheme for calculating the initial time derivatives,
  // but the rhs of the standard element routine is used as starting point for this special system
  // of equations. Therefore, the rhs vector has to be scaled correctly.
  SetElementTimeParameter(true);

  // deactivate turbulence settings
  SetElementTurbulenceParameters(true);

  // call core algorithm
  ScaTraTimIntImpl::CalcInitialTimeDerivative();

  // and finally undo our temporary settings
  SetElementGeneralParameters(false);
  SetElementTimeParameter(false);
  SetElementTurbulenceParameters(false);

  return;
}


/*--------------------------------------------------------------------*
 | set state on micro scale in multi-scale simulations     fang 11/15 |
 *--------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetState(
    Teuchos::RCP<Epetra_Vector> phin,   //!< micro-scale state vector at old time step
    Teuchos::RCP<Epetra_Vector> phinp,  //!< micro-scale state vector at new time step
    Teuchos::RCP<Epetra_Vector>
        phidtn,  //!< time derivative of micro-scale state vector at old time step
    Teuchos::RCP<Epetra_Vector>
        phidtnp,  //!< time derivative of micro-scale state vector at new time step
    Teuchos::RCP<Epetra_Vector> hist,               //!< micro-scale history vector
    Teuchos::RCP<IO::DiscretizationWriter> output,  //!< micro-scale discretization writer
    const std::vector<double>&
        phinp_macro,   //!< values of state variables at macro-scale Gauss point
    const int step,    //!< time step
    const double time  //!< time
)
{
  phin_ = phin;
  phinp_ = phinp;
  phidtn_ = phidtn;
  phidtnp_ = phidtnp;
  hist_ = hist;
  output_ = output;
  phinp_macro_ = phinp_macro;
  dq_dphi_.resize(phinp_macro_.size(), 0.);
  step_ = step;
  time_ = time;

  return;
}  // SCATRA::TimIntOneStepTheta::SetState


/*--------------------------------------------------------------------*
 | clear state on micro scale in multi-scale simulations   fang 11/15 |
 *--------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ClearState()
{
  phin_ = Teuchos::null;
  phinp_ = Teuchos::null;
  phidtn_ = Teuchos::null;
  phidtnp_ = Teuchos::null;
  hist_ = Teuchos::null;
  output_ = Teuchos::null;
  phinp_macro_.clear();
  dq_dphi_.clear();
  step_ = -1;
  time_ = 0.;

  return;
}  // SCATRA::TimIntOneStepTheta::ClearState
