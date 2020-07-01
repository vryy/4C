/*----------------------------------------------------------------------*/
/*! \file
\brief BDF2 time-integration scheme

\level 1


*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_meshtying_strategy_base.H"
#include "turbulence_hit_scalar_forcing.H"

#include "../drt_io/io.H"

#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "scatra_timint_bdf2.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntBDF2::TimIntBDF2(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, params, extraparams, output),
      theta_(1.0),
      phinm_(Teuchos::null),
      fsphinp_(Teuchos::null)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                         rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::Setup()
{
  // initialize base class
  ScaTraTimIntImpl::Setup();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // state vector for solution at time t_{n-1}
  phinm_ = LINALG::CreateVector(*dofrowmap, true);

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

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   gjb 08/08 |
*----------------------------------------------------------------------*/
SCATRA::TimIntBDF2::~TimIntBDF2() { return; }


/*----------------------------------------------------------------------*
 |  set time parameter for element evaluation (usual call)   ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::SetElementTimeParameter(bool forcedincrementalsolver) const
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
void SCATRA::TimIntBDF2::SetTimeForNeumannEvaluation(Teuchos::ParameterList& params)
{
  params.set("total time", time_);
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::SetOldPartOfRighthandside()
{
  // call base class routine
  ScaTraTimIntImpl::SetOldPartOfRighthandside();

  /*
  BDF2: for variable time step:

                 hist_ = (1+omega)^2/(1+ 2*omega) * phin_
                           - omega^2/(1+ 2*omega) * phinm_

  BDF2: for constant time step:

                 hist_ = 4/3 phin_ - 1/3 phinm_
  */
  if (step_ > 1)
  {
    double fact1 = 4.0 / 3.0;
    double fact2 = -1.0 / 3.0;
    hist_->Update(fact1, *phin_, fact2, *phinm_, 0.0);

    // for BDF2 theta is set to 2/3 for constant time-step length dt
    theta_ = 2.0 / 3.0;
  }
  else
  {
    // for start-up of BDF2 we do one step with backward Euler
    hist_->Update(1.0, *phin_, 0.0);

    // backward Euler => use theta=1.0
    theta_ = 1.0;
  }

  return;
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                         gjb 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::ExplicitPredictor() const
{
  // call base class routine
  ScaTraTimIntImpl::ExplicitPredictor();

  if (step_ > 1) phinp_->Update(-1.0, *phinm_, 2.0);
  // for step == 1 phinp_ is already correctly initialized with the
  // initial field phin_

  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                             |
 | scaled with a factor resulting from time discretization     vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::AddNeumannToResidual()
{
  residual_->Update(theta_ * dta_, *neumann_loads_, 1.0);
  return;
}


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::AVM3Separation()
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
void SCATRA::TimIntBDF2::DynamicComputationOfCs()
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
 | dynamic Vreman model                                krank      09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::DynamicComputationOfCv()
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
void SCATRA::TimIntBDF2::AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver)
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
void SCATRA::TimIntBDF2::ComputeTimeDerivative()
{
  // call base class routine
  ScaTraTimIntImpl::ComputeTimeDerivative();

  if (step_ == 1)
  {
    // time derivative of phi for first time step:
    // phidt(n+1) = (phi(n+1)-phi(n))/dt
    const double fact = 1.0 / dta_;
    phidtnp_->Update(fact, *phinp_, -fact, *hist_, 0.0);
  }
  else
  {
    // time derivative of phi:
    // phidt(n+1) = ((3/2)*phi(n+1)-2*phi(n)+(1/2)*phi(n-1))/dt
    const double fact = 3.0 / (2.0 * dta_);
    phidtnp_->Update(fact, *phinp_, -fact, *hist_, 0.0);
  }

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
void SCATRA::TimIntBDF2::Update(const int num)
{
  // call base class routine
  ScaTraTimIntImpl::Update(num);

  // compute flux vector field for later output BEFORE time shift of results
  // is performed below !!
  if (calcflux_domain_ != INPAR::SCATRA::flux_none or
      calcflux_boundary_ != INPAR::SCATRA::flux_none)
  {
    if (DoOutput() or DoBoundaryFluxStatistics()) CalcFlux(true, num);
  }

  // solution of this step becomes most recent solution of the last step
  phinm_->Update(1.0, *phin_, 0.0);
  phin_->Update(1.0, *phinp_, 0.0);

  // call time update of forcing routine
  if (homisoturb_forcing_ != Teuchos::null) homisoturb_forcing_->TimeUpdateForcing();

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::OutputRestart() const
{
  // call base class routine
  ScaTraTimIntImpl::OutputRestart();

  // additional state vectors that are needed for BDF2 restart
  output_->WriteVector("phin", phin_);
  output_->WriteVector("phinm", phinm_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntBDF2::ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input)
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

  // read state vectors that are needed for BDF2 restart
  reader->ReadVector(phinp_, "phinp");
  reader->ReadVector(phin_, "phin");
  reader->ReadVector(phinm_, "phinm");

  ReadRestartProblemSpecific(step, *reader);

  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no or
      turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    AVM3Preparation();

  return;
}
