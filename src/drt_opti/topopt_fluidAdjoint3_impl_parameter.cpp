/*---------------------------------------------------------------------*/
/*! \file

\brief general element parameter for fluid adjoint equations for topology optimization

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#include "topopt_fluidAdjoint3_impl_parameter.H"
#include "../drt_lib/drt_dserror.H"


//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::FluidAdjoint3ImplParameter>
    DRT::ELEMENTS::FluidAdjoint3ImplParameter::instance_;


//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::FluidAdjoint3ImplParameter>
DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance()
{
  if (instance_ == Teuchos::null) instance_ = Teuchos::rcp(new FluidAdjoint3ImplParameter());
  return instance_;
}


//----------------------------------------------------------------------*/
// private constructor of FluidImplParameter
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidAdjoint3ImplParameter::FluidAdjoint3ImplParameter()
    : set_general_adjoint_parameter_(false),
      dissipation_(INPAR::TOPOPT::obj_diss_no),
      pressure_drop_(false),
      dissipation_fac_(0.0),
      pressure_drop_fac_(0.0),
      adjoint_type_(INPAR::TOPOPT::discrete_adjoint),
      is_stationary_(false),
      is_init_instat_step_(false),
      is_inconsistent_(false),
      pspg_(true),
      supg_(true),
      graddiv_(true),
      whichtau_(INPAR::FLUID::tau_not_defined),
      mat_gp_(false),  // standard evaluation of the material at the element center
      tau_gp_(false),  // standard evaluation of tau at the element center
      timealgo_(INPAR::FLUID::timeint_one_step_theta),
      time_(-1.0),
      dt_(0.0),
      theta_obj_(0.0),
      theta_(0.0),
      omtheta_(0.0),
      timefac_(0.0),
      timefacrhs_(0.0)
{
}


//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidAdjoint3ImplParameter::SetElementGeneralAdjointParameter(
    Teuchos::ParameterList& params)
{
  if (set_general_adjoint_parameter_ == false) set_general_adjoint_parameter_ = true;

  //----------------------------------------------------------------------
  // get flags to switch on/off different fluid formulations
  //----------------------------------------------------------------------

  // set flag, time integration scheme
  timealgo_ = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "TimeIntegrationScheme");

  adjoint_type_ = params.get<INPAR::TOPOPT::AdjointType>("adjoint type");

  // set time integration scheme-specific element parameters
  if (timealgo_ == INPAR::FLUID::timeint_stationary)
    is_stationary_ = true;
  else
    is_stationary_ = false;

  // set if objective contains dissipation and according factor
  dissipation_ = params.get<INPAR::TOPOPT::ObjectiveDissipation>("dissipation");
  if (dissipation_ != INPAR::TOPOPT::obj_diss_no)
    dissipation_fac_ = params.get<double>("dissipation_fac");
  else
    dissipation_fac_ = 0.0;

  // set if objective contains pressure drop and according factor
  pressure_drop_ = params.get<bool>("presDrop");
  if (pressure_drop_)
    pressure_drop_fac_ = params.get<double>("presDropFac");
  else
    pressure_drop_fac_ = 0.0;


  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------
  Teuchos::ParameterList& stablist = params.sublist("RESIDUAL-BASED STABILIZATION");

  // no safety check necessary since all options are used
  pspg_ = DRT::INPUT::IntegralValue<int>(stablist, "PSPG");
  supg_ = DRT::INPUT::IntegralValue<int>(stablist, "SUPG");
  graddiv_ = DRT::INPUT::IntegralValue<int>(stablist, "GRAD_DIV");

  //-------------------------------
  // get tau definition
  //-------------------------------
  whichtau_ = DRT::INPUT::IntegralValue<INPAR::FLUID::TauType>(stablist, "DEFINITION_TAU");
  // check if tau can be handled
  if (not(whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins or
          whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt or
          whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen or
          whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt or
          whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins_scaled or
          whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt or
          whichtau_ == INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall or
          whichtau_ == INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt or
          whichtau_ == INPAR::FLUID::tau_shakib_hughes_codina or
          whichtau_ == INPAR::FLUID::tau_shakib_hughes_codina_wo_dt or
          whichtau_ == INPAR::FLUID::tau_codina or whichtau_ == INPAR::FLUID::tau_codina_wo_dt or
          whichtau_ == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina or
          whichtau_ == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt))
    dserror("Definition of Tau cannot be handled by the element");

  // set correct stationary definition of stabilization parameter automatically
  if (is_stationary_)
  {
    if (whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins)
      whichtau_ = INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt;
    else if (whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen)
      whichtau_ = INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt;
    else if (whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins_scaled)
      whichtau_ = INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt;
    else if (whichtau_ == INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall)
      whichtau_ = INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt;
    else if (whichtau_ == INPAR::FLUID::tau_shakib_hughes_codina)
      whichtau_ = INPAR::FLUID::tau_shakib_hughes_codina_wo_dt;
    else if (whichtau_ == INPAR::FLUID::tau_codina)
      whichtau_ = INPAR::FLUID::tau_codina_wo_dt;
    else if (whichtau_ == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina)
      whichtau_ = INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt;
  }

  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (stablist.get<std::string>("STABTYPE") == "inconsistent") is_inconsistent_ = true;


  // set flags for potential evaluation of tau and material law at int. point
  // default value: evaluation at element center
  const std::string tauloc = stablist.get<std::string>("EVALUATION_TAU");
  if (tauloc == "integration_point")
    tau_gp_ = true;
  else
    tau_gp_ = false;
  const std::string matloc = stablist.get<std::string>("EVALUATION_MAT");
  if (matloc == "integration_point")
    mat_gp_ = true;
  else
    mat_gp_ = false;

  //-------------------------------
  // get material parameter
  //-------------------------------
  dens_ = params.get<double>("density");
  visc_ = params.get<double>("viscosity");

  topopt_params_[0] = params.get<double>("MIN_PORO");
  topopt_params_[1] = params.get<double>("MAX_PORO");
  topopt_params_[2] = params.get<double>("SMEAR_FAC");

  // set flag for test cases
  testcase_ = params.get<INPAR::TOPOPT::AdjointCase>("special test case");
}

void DRT::ELEMENTS::FluidAdjoint3ImplParameter::SetElementAdjointTimeParameter(
    Teuchos::ParameterList& params)
{
  // second check: timealgo
  // work around to use SetTimeParameter in GenaAlpha (Neumann BC)
  if (set_general_adjoint_parameter_ != true)
    dserror("General adjoint parameter are not set yet!!");

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------

  // get current time: n+alpha_F for generalized-alpha scheme, n+1 otherwise
  time_ = params.get<double>("total time", -1.0);

  // get flag if step 0 in instationary
  is_init_instat_step_ = params.get<bool>("initial instationary step", false);

  // get time-step length and time-integration parameters
  if (not is_stationary_)
  {
    dt_ = params.get<double>("dt");

    theta_obj_ = params.get<double>("theta_obj");

    theta_ = params.get<double>("theta");
    omtheta_ = params.get<double>("omtheta");
  }
  else  // is_stationary == true
  {
    theta_obj_ = 1.0;

    dt_ = theta_ = 1.0;
    omtheta_ = 0.0;
  }

  timefac_ = theta_ * dt_;
  timefacrhs_ = omtheta_ * dt_;


  double TOL = 1.0e-14;
  if (dt_ < 0.0 or time_ < -TOL or theta_obj_ < 0.0 or theta_obj_ > 1.0 or theta_ < 0.0 or
      omtheta_ < 0.0)
  {
    std::cout << "dt_: " << dt_ << std::endl;
    std::cout << "time_ " << time_ << std::endl;
    std::cout << "theta_obj_ " << theta_obj_ << std::endl;
    std::cout << "theta_ " << theta_ << std::endl;
    std::cout << "omtheta_ " << omtheta_ << std::endl;
    dserror("Negative (or no) time-integration parameter or time-step length supplied");
  }
}



//----------------------------------------------------------------------*/
// print fluid parameter to screen (AE 01-11)
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidAdjoint3ImplParameter::PrintAdjointParameter() const
{
  std::cout << std::endl
            << "|-----------------------------------------------------------------------------"
            << std::endl;
  std::cout << "|  General Fluid Adjoint parameter: " << std::endl;
  std::cout << "|-----------------------------------------------------------------------------"
            << std::endl;
  //! flag if continous or discrete adjoints are computed
  std::cout << "|    type of adjoint equations is:    " << adjoint_type_ << std::endl;
  //! flag SetGeneralParameter was called
  std::cout << "|    method SetElmentGeneralAdjointParameter was called:    "
            << set_general_adjoint_parameter_ << std::endl;
  //! flag to (de)activate stationary formulation
  std::cout << "|    steady state:    " << is_stationary_ << std::endl;
  //! flag to (de)activate special treatment of step 0 in instationary computation
  std::cout << "|    initial instationary step:    " << is_init_instat_step_ << std::endl;
  //! flag to (de)activate second derivatives
  std::cout << "|    use inconsistent:    " << is_inconsistent_ << std::endl;
  //! Flag to (de)activate PSPG stabilization
  std::cout << "|    PSPG:    " << pspg_ << std::endl;
  //! Flag to (de)activate SUPG stabilization
  std::cout << "|    SUPG:    " << supg_ << std::endl;
  //! Flag to (de)activate least-squares stabilization of continuity equation
  std::cout << "|    Grad-Div-Stab:    " << graddiv_ << std::endl;
  //! Flag to define tau
  std::cout << "|    Definition of stabilization parameter:    " << whichtau_ << std::endl;
  //! flag for material evaluation at Gaussian integration points
  std::cout << "|    material evaluation at Gaussian integration points:    " << mat_gp_
            << std::endl;
  //! flag for stabilization parameter evaluation at Gaussian integration points
  std::cout << "|    stabilization parameter evaluation at Gaussian integration points:    "
            << tau_gp_ << std::endl;
  //! value of physical density
  std::cout << "|    physical density:    " << dens_ << std::endl;
  //! value of physical viscosity
  std::cout << "|    physical viscosity:    " << visc_ << std::endl;
  //! matrix with values for computation of porosity with respect to topopt density
  std::cout << "|    topology optimization porosity parameter: " << topopt_params_[0] << " "
            << topopt_params_[1] << " " << topopt_params_[2] << std::endl;
  /// enumeration of testcase
  std::cout << "|    academical testcase:    " << testcase_ << std::endl;
  std::cout << "|---------------------------------------------------------------------------"
            << std::endl;

  std::cout << std::endl
            << "|---------------------------------------------------------------------------"
            << std::endl;
  std::cout << "|  Time parameter: " << std::endl;
  std::cout << "|---------------------------------------------------------------------------"
            << std::endl;
  //! time algorithm
  std::cout << "|    time algorithm:    " << timealgo_ << std::endl;
  //! actual time to evaluate the body BC
  std::cout << "|    time:    " << time_ << std::endl;
  //! time-step length
  std::cout << "|    time step:    " << dt_ << std::endl;
  //! factor for right-hand side due to one-step-theta time-integration scheme of objective function
  std::cout << "|    theta_obj:    " << theta_obj_ << std::endl;
  //! factor for left-hand side due to one-step-theta time-integration scheme
  std::cout << "|    theta:    " << theta_ << std::endl;
  //! factor for right-hand side due to one-step-theta time-integration scheme
  std::cout << "|    (1-theta):    " << omtheta_ << std::endl;
  //! timefac = dt_ * theta_
  std::cout << "|    time factor:    " << timefac_ << std::endl;
  //! timefacrhs = dt_ * (1-theta_)
  std::cout << "|    time factor rhs:    " << timefacrhs_ << std::endl;
  std::cout << "|---------------------------------------------------------------------------"
            << std::endl;
}
