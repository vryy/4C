/*---------------------------------------------------------------------*/
/*!

\brief Evaluation of element parameter

\level 3

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include "topopt_optimizer_ele_parameter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"



//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::TopOptParam> DRT::ELEMENTS::TopOptParam::instance_;


//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::TopOptParam> DRT::ELEMENTS::TopOptParam::Instance()
{
  if (instance_ == Teuchos::null) instance_ = Teuchos::rcp(new TopOptParam());
  return instance_;
}


//----------------------------------------------------------------------*/
// private constructor of FluidImplParameter
//----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptParam::TopOptParam()
    : dens_(-1.0),
      visc_(-1.0),
      min_poro_(-1.0),
      max_poro_(-1.0),
      smear_fac_(-1.0),
      dissipation_(INPAR::TOPOPT::obj_diss_no),
      pressure_drop_(false),
      dissipation_fac_(0.0),
      pressure_drop_fac_(0.0),
      is_stationary_(false),
      timealgo_(INPAR::FLUID::timeint_one_step_theta),
      supg_(false),
      pspg_(false),
      whichtau_(INPAR::FLUID::tau_not_defined),
      dt_(-1.0),
      max_timesteps_(-1),
      theta_obj_(-1.0),
      theta_(-1.0),
      vol_bd_(-1.0),
      dens_type_(INPAR::TOPOPT::dens_undefined),
      opti_case_(INPAR::TOPOPT::optitest_no)
{
}


//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOptParam::SetGeneralOptimizationParameter(Teuchos::ParameterList& params)
{
  //----------------------------------------------------------------------
  // get flags to switch on/off different fluid formulations
  //----------------------------------------------------------------------

  // flow material parameter
  dens_ = params.get<double>("density");
  visc_ = params.get<double>("viscosity");

  // optimization material parameter
  min_poro_ = params.get<double>("MIN_PORO");
  max_poro_ = params.get<double>("MAX_PORO");
  smear_fac_ = params.get<double>("SMEAR_FAC");

  // check whether there is zero or negative (physical) viscosity
  // (expect for permeable fluid)
  if (visc_ < EPS15) dserror("zero or negative (physical) diffusivity");

  // set flag, time integration scheme
  timealgo_ = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");

  // set time integration scheme-specific element parameters
  if (timealgo_ == INPAR::FLUID::timeint_stationary)
    is_stationary_ = true;
  else
    is_stationary_ = false;

  supg_ = DRT::INPUT::IntegralValue<int>(params.sublist("RESIDUAL-BASED STABILIZATION"), "SUPG");
  pspg_ = DRT::INPUT::IntegralValue<int>(params.sublist("RESIDUAL-BASED STABILIZATION"), "PSPG");

  //-------------------------------
  // get tau definition
  //-------------------------------
  whichtau_ = DRT::INPUT::IntegralValue<INPAR::FLUID::TauType>(
      params.sublist("RESIDUAL-BASED STABILIZATION"), "DEFINITION_TAU");
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

  // set if objective contains dissipation and according factor
  dissipation_ = params.get<INPAR::TOPOPT::ObjectiveDissipation>("dissipation");
  if (dissipation_ != INPAR::TOPOPT::obj_diss_no)
    dissipation_fac_ = params.get<double>("dissipation_fac");
  else
    dissipation_fac_ = 0.0;

  // set if objective contains pressure drop and according factor
  pressure_drop_ = params.get<bool>("pres_drop");
  if (pressure_drop_)
    pressure_drop_fac_ = params.get<double>("pres_drop_fac");
  else
    pressure_drop_fac_ = 0.0;

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------
  if (not is_stationary_)
  {
    theta_obj_ = params.get<double>("theta_obj");

    dt_ = params.get<double>("dt");
    max_timesteps_ = params.get<int>("maxtimesteps");
    theta_ = params.get<double>("theta");
  }
  else
  {
    theta_obj_ = 1.0;

    dt_ = theta_ = 1.0;
    max_timesteps_ = 1;
  }

  vol_bd_ = params.get<double>("vol_bd");

  dens_type_ = params.get<INPAR::TOPOPT::DensityField>("dens_type");

  opti_case_ = params.get<INPAR::TOPOPT::OptiCase>("opti_case");
}



//----------------------------------------------------------------------*
//  update general parameters                          winklmaier 12/14 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOptParam::UpdateGeneralOptimizationParameter(Teuchos::ParameterList& params)
{
  smear_fac_ = params.get<double>("SMEAR_FAC");
}



//----------------------------------------------------------------------*/
// print fluid parameter to screen (AE 01-11)
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOptParam::PrintAdjointParameter() const
{
  std::cout << std::endl
            << "|-----------------------------------------------------------------------------"
            << std::endl;
  std::cout << "|  Material parameter: " << std::endl;
  std::cout << "|-----------------------------------------------------------------------------"
            << std::endl;
  // boolean is true if objective contains dissipation
  std::cout << "|    physical density    " << dens_ << std::endl;
  // boolean is true if objective contains pressure drop
  std::cout << "|    physical viscosity    " << visc_ << std::endl;
  // minimal inverse modelling porosity
  std::cout << "|    minimal pseudo-porosity:    " << min_poro_ << std::endl;
  // maximal inverse modelling porosity
  std::cout << "|    maximal pseudo-porosity:    " << max_poro_ << std::endl;
  // smearing factor between density and porosity
  std::cout << "|    smearing factor:    " << smear_fac_ << std::endl;

  std::cout << "|  General optimization parameter: " << std::endl;
  std::cout << "|-----------------------------------------------------------------------------"
            << std::endl;
  // boolean is true if objective contains dissipation
  std::cout << "|    objective dissipation on?    " << dissipation_ << std::endl;
  // boolean is true if objective contains pressure drop
  std::cout << "|    objective pressure drop on?    " << pressure_drop_ << std::endl;
  // objective's dissipation factor
  std::cout << "|    objective dissipation factor:    " << dissipation_fac_ << std::endl;
  // objective's pressure drop factor
  std::cout << "|    objective pressure drop factor:    " << pressure_drop_fac_ << std::endl;
  // type of optimization field
  std::cout << "|    type of optimization field:    " << dens_type_ << std::endl;
  // test case scenario
  std::cout << "|    optimization test case number:    " << opti_case_ << std::endl;

  std::cout << std::endl
            << "|---------------------------------------------------------------------------"
            << std::endl;
  std::cout << "|  Flow parameter: " << std::endl;
  std::cout << "|---------------------------------------------------------------------------"
            << std::endl;
  //! flag to (de)activate stationary formulation
  std::cout << "|    steady state:    " << is_stationary_ << std::endl;
  //! time algorithm
  std::cout << "|    time algorithm:    " << timealgo_ << std::endl;
  //! Flag to define tau
  std::cout << "|    Definition of stabilization parameter:    " << whichtau_ << std::endl;
  //! time-step length
  std::cout << "|    time step:    " << dt_ << std::endl;
  /// maximal number of time steps
  std::cout << "|    maximal number of time steps:     " << max_timesteps_ << std::endl;
  /// theta for objective integration
  std::cout << "|    theta_obj:     " << theta_obj_ << std::endl;
  /// theta
  std::cout << "|    theta:     " << theta_ << std::endl;
  std::cout << "|---------------------------------------------------------------------------"
            << std::endl;
}  /// @name objective parameters
