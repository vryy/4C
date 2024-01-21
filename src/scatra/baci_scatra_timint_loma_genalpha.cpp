/*----------------------------------------------------------------------*/
/*! \file

\brief Generalized-alpha time-integration scheme with extensions for
       loma problems

\level 2

*/
/*----------------------------------------------------------------------*/

#include "baci_scatra_timint_loma_genalpha.H"

#include "baci_fluid_turbulence_dyn_smag.H"
#include "baci_fluid_turbulence_dyn_vreman.H"
#include "baci_global_data.H"
#include "baci_io.H"
#include "baci_lib_utils_parameter_list.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_scatra_ele_action.H"
#include "baci_scatra_turbulence_hit_scalar_forcing.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       vg 11/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntLomaGenAlpha::TimIntLomaGenAlpha(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntLoma(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntGenAlpha(actdis, solver, sctratimintparams, extraparams, output),
      thermpressaf_(0.0),
      thermpressam_(0.0),
      thermpressdtaf_(0.0),
      thermpressdtam_(0.0)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::Init();
  ScaTraTimIntLoma::Init();

  return;
}

/*----------------------------------------------------------------------*
 |  setup time integration                                  rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::Setup()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::Setup();
  ScaTraTimIntLoma::Setup();

  return;
}



/*----------------------------------------------------------------------*
 | predict thermodynamic pressure and time derivative          vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::PredictThermPressure()
{
  // same-thermodynamic-pressure predictor (not required to be performed,
  // since we just updated the thermodynamic pressure, and thus,
  // thermpressnp_ = thermpressn_)
  // prediction of time derivative:
  double fact = (gamma_ - 1.0) / gamma_;
  thermpressdtnp_ = fact * thermpressdtn_;

  // same-thermodynamic-pressure-derivative predictor (currrently not used)
  // thermpressnp_ += dta_*thermpressdtn_;
  // prediction of time derivative not required (would also not be required
  // to be performed, since we just updated the time derivatives of density,
  // and thus, thermpressdtnp_ = thermpressdtn_)

  return;
}


/*----------------------------------------------------------------------*
 | compute values of therm. pressure at interm. time steps     vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::ComputeThermPressureIntermediateValues()
{
  // thermodynamic pressure at n+alpha_F and n+alpha_M for low-Mach-number case
  // -> required for evaluation of equation of state
  thermpressaf_ = alphaF_ * thermpressnp_ + (1.0 - alphaF_) * thermpressn_;
  thermpressam_ = alphaM_ * thermpressnp_ + (1.0 - alphaM_) * thermpressn_;

  // time derivative of thermodyn. press. at n+alpha_F for low-Mach-number case
  // -> required as right-hand-side contribution to temperature equation,
  // hence, evaluated at n+alpha_F
  thermpressdtaf_ = alphaF_ * thermpressdtnp_ + (1.0 - alphaF_) * thermpressdtn_;

  // time derivative of thermodyn. press. at n+alpha_M for low-Mach-number case
  // -> required for transfer to flow solver and use in continuity equation
  thermpressdtam_ = alphaM_ * thermpressdtnp_ + (1.0 - alphaM_) * thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | compute thermodynamic pressure for low-Mach-number flow     vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::ComputeThermPressure()
{
  // compute temperature at n+alpha_F
  phiaf_->Update(alphaF_, *phinp_, (1.0 - alphaF_), *phin_, 0.0);

  // define element parameter list
  Teuchos::ParameterList eleparams;

  // DO THIS BEFORE PHINP IS SET (ClearState() is called internally!!!!)
  // compute flux approximation and add it to the parameter list
  AddFluxApproxToParameterList(eleparams);

  // set scalar values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp", phiaf_);

  // set action for elements
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_domain_and_bodyforce, eleparams);
  SetElementTimeParameter();

  // variables for integrals of domain and bodyforce
  Teuchos::RCP<CORE::LINALG::SerialDenseVector> scalars =
      Teuchos::rcp(new CORE::LINALG::SerialDenseVector(2));

  // evaluate domain and bodyforce integral
  discret_->EvaluateScalars(eleparams, scalars);

  // get global integral values
  double pardomint = (*scalars)[0];
  double parbofint = (*scalars)[1];

  // set action for elements
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
      "action", SCATRA::BoundaryAction::calc_loma_therm_press, eleparams);

  // variables for integrals of normal velocity and diffusive flux
  double normvelint = 0.0;
  double normdifffluxint = 0.0;
  eleparams.set("normal velocity integral", normvelint);
  eleparams.set("normal diffusive flux integral", normdifffluxint);

  // evaluate velocity-divergence and diffusive (minus sign!) flux on boundaries
  // We may use the flux-calculation condition for calculation of fluxes for
  // thermodynamic pressure, since it is usually at the same boundary.
  std::vector<std::string> condnames;
  condnames.push_back("ScaTraFluxCalc");
  for (unsigned int i = 0; i < condnames.size(); i++)
  {
    discret_->EvaluateCondition(eleparams, Teuchos::null, Teuchos::null, Teuchos::null,
        Teuchos::null, Teuchos::null, condnames[i]);
  }

  // get integral values on this proc
  normvelint = eleparams.get<double>("normal velocity integral");
  normdifffluxint = eleparams.get<double>("normal diffusive flux integral");

  // get integral values in parallel case
  double parnormvelint = 0.0;
  double parnormdifffluxint = 0.0;
  discret_->Comm().SumAll(&normvelint, &parnormvelint, 1);
  discret_->Comm().SumAll(&normdifffluxint, &parnormdifffluxint, 1);

  // clean up
  discret_->ClearState();

  // compute thermodynamic pressure (with specific heat ratio fixed to be 1.4)
  const double shr = 1.4;
  const double divt = shr * parnormvelint / pardomint;
  const double lhs = alphaF_ * genalphafac_ * dta_ * divt;
  const double rhs =
      genalphafac_ * dta_ * (shr - 1.0) * (-parnormdifffluxint + parbofint) / pardomint;
  const double hist = thermpressn_ - (1.0 - alphaF_) * genalphafac_ * dta_ * divt * thermpressn_ +
                      (1.0 - genalphafac_) * dta_ * thermpressdtn_;
  thermpressnp_ = (rhs + hist) / (1.0 + lhs);

  // print out thermodynamic pressure
  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
    std::cout << "Data output for instationary thermodynamic pressure:" << std::endl;
    std::cout << "Velocity in-/outflow at indicated boundary: " << parnormvelint << std::endl;
    std::cout << "Diffusive flux at indicated boundary: " << parnormdifffluxint << std::endl;
    std::cout << "Thermodynamic pressure: " << thermpressnp_ << std::endl;
    std::cout << "+--------------------------------------------------------------------------------"
                 "------------+"
              << std::endl;
  }

  // compute time derivative of thermodynamic pressure at time step n+1
  ComputeThermPressureTimeDerivative();

  // compute values at intermediate time steps
  ComputeThermPressureIntermediateValues();

  return;
}


/*----------------------------------------------------------------------*
 | compute time derivative of thermodynamic pressure           vg 09/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::ComputeThermPressureTimeDerivative()
{
  // time derivative of thermodynamic pressure:
  // tpdt(n+1) = (tp(n+1)-tp(n)) / (gamma*dt) + (1-(1/gamma))*tpdt(n)
  const double fact1 = 1.0 / (gamma_ * dta_);
  const double fact2 = 1.0 - (1.0 / gamma_);
  thermpressdtnp_ = fact1 * (thermpressnp_ - thermpressn_) + fact2 * thermpressdtn_;

  return;
}


/*----------------------------------------------------------------------*
 | update thermodynamic pressure at n for low-Mach-number flow vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::UpdateThermPressure()
{
  thermpressn_ = thermpressnp_;
  thermpressdtn_ = thermpressdtnp_;

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::WriteRestart() const
{
  // write standard fields first
  TimIntGenAlpha::WriteRestart();

  // write additional restart data for loma
  // required for restart of closed systems

  // thermodynamic pressure at time n+1
  output_->WriteDouble("thermpressnp", thermpressnp_);
  // thermodynamic pressure at time n
  output_->WriteDouble("thermpressn", thermpressn_);
  // thermodynamic pressure at time n+alpha_f
  output_->WriteDouble("thermpressaf", thermpressaf_);
  // thermodynamic pressure at time n+alpha_m
  output_->WriteDouble("thermpressam", thermpressam_);
  // time derivative of thermodynamic pressure at time n+1
  output_->WriteDouble("thermpressdtnp", thermpressdtnp_);
  // time derivative of thermodynamic pressure at time n
  output_->WriteDouble("thermpressdtn", thermpressdtn_);
  // time derivative of thermodynamic pressure at time n+alpha_f
  output_->WriteDouble("thermpressdtaf", thermpressdtaf_);
  // time derivative of thermodynamic pressure at time n+alpha_m
  output_->WriteDouble("thermpressdtam", thermpressdtam_);
  // as well as initial mass
  output_->WriteDouble("initialmass", initialmass_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                             vg 11/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input)
{
  // do standard output
  TimIntGenAlpha::ReadRestart(step, input);

  // restart data of loma problems
  // required for restart of closed systems

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // thermodynamic pressure at time n+1
  thermpressnp_ = reader->ReadDouble("thermpressnp");
  // thermodynamic pressure at time n
  thermpressn_ = reader->ReadDouble("thermpressn");
  // thermodynamic pressure at time n+alpha_f
  thermpressaf_ = reader->ReadDouble("thermpressaf");
  // thermodynamic pressure at time n+alpha_m
  thermpressam_ = reader->ReadDouble("thermpressam");
  // time derivative of thermodynamic pressure at time n+1
  thermpressdtnp_ = reader->ReadDouble("thermpressdtnp");
  // time derivative of thermodynamic pressure at time n
  thermpressdtn_ = reader->ReadDouble("thermpressdtn");
  // time derivative of thermodynamic pressure at time n+alpha_f
  thermpressdtaf_ = reader->ReadDouble("thermpressdtaf");
  // time derivative of thermodynamic pressure at time n+alpha_m
  thermpressdtam_ = reader->ReadDouble("thermpressdtam");
  // as well as initial mass
  initialmass_ = reader->ReadDouble("initialmass");

  return;
}


/*----------------------------------------------------------------------*
 | dynamic Smagorinsky model                           rasthofer  08/12 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::DynamicComputationOfCs()
{
  if (turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
  {
    // perform filtering and computation of Prt
    // compute averaged values for LkMk and MkMk
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    DynSmag_->ApplyFilterForDynamicComputationOfPrt(
        phiaf_, thermpressaf_, dirichtoggle, *extraparams_, NdsVel());
  }

  return;
}


/*----------------------------------------------------------------------*
 | dynamic Vreman model                                krank  09/13     |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::DynamicComputationOfCv()
{
  if (turbmodel_ == INPAR::FLUID::dynamic_vreman)
  {
    const Teuchos::RCP<const Epetra_Vector> dirichtoggle = DirichletToggle();
    Vrem_->ApplyFilterForDynamicComputationOfDt(
        phiaf_, thermpressaf_, dirichtoggle, *extraparams_, NdsVel());
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 | add thermodynamic pressure to parameter list for element evaluation rasthofer 12/13 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::TimIntLomaGenAlpha::AddThermPressToParameterList(
    Teuchos::ParameterList& params  //!< parameter list
)
{
  params.set("thermodynamic pressure", thermpressaf_);
  params.set("thermodynamic pressure at n+alpha_M", thermpressam_);
  params.set("time derivative of thermodynamic pressure", thermpressdtaf_);
  discret_->SetState("phiam", phiam_);
  return;
}

BACI_NAMESPACE_CLOSE
