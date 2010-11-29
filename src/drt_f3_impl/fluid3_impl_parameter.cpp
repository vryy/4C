/*----------------------------------------------------------------------*/
/*!
\file fluid3_impl_parameter.cpp

\brief Evaluation of general fluid parameter

Fluid3ImplParameter::SetParameter(Teuchos::ParameterList& params)
set all general fluid parameter once for all elements.

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include <string>

#include "fluid3_impl_parameter.H"
#include "../drt_lib/drt_dserror.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3ImplParameter* DRT::ELEMENTS::Fluid3ImplParameter::instance_;


//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3ImplParameter* DRT::ELEMENTS::Fluid3ImplParameter::Instance()
{
  if (instance_==NULL)
    instance_ = new Fluid3ImplParameter();
  return instance_;
}


//----------------------------------------------------------------------*/
// private constructor of Fluid3ImplParameter
//----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3ImplParameter::Fluid3ImplParameter()
  :
  is_genalpha_(false),
  is_conservative_(false),
  is_stationary_(false),
  is_newton_(false),
  is_inconsistent_(false),
  physicaltype_(INPAR::FLUID::incompressible),
  tds_(INPAR::FLUID::subscales_none),
  transient_(INPAR::FLUID::inertia_stab_drop),
  pspg_(INPAR::FLUID::pstab_use_pspg),
  supg_(INPAR::FLUID::convective_stab_supg),
  vstab_(INPAR::FLUID::viscous_stab_none),
  cstab_(INPAR::FLUID::continuity_stab_yes),
  cross_(INPAR::FLUID::cross_stress_stab_none),
  reynolds_(INPAR::FLUID::reynolds_stress_stab_none),
  whichtau_(INPAR::FLUID::tau_not_defined),
  fssgv_(INPAR::FLUID::no_fssgv),
  turb_mod_action_(INPAR::FLUID::no_model),
  mat_gp_(false),     // standard evaluation of the material at the element center
  tau_gp_(false),     // standard evaluation of tau at the element center
  timefac_(0.0),
  Cs_(0.0),
  l_tau_(0.0)
{
}


//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3ImplParameter::SetParameter( Teuchos::ParameterList& params )
{
//----------------------------------------------------------------------
// get flags to switch on/off different fluid formulations
//----------------------------------------------------------------------

  // set flag, if it is a stationary fluid system
  is_stationary_ = params.get<bool>("is stationary");

  // check whether we have a generalized-alpha time-integration scheme
  is_genalpha_ = params.get<bool>("using generalized-alpha time integration");

  // set flag for type of linearization (fixed-point-like or Newton)
  //std::string newtonstr   = params.get<std::string>("Linearisation");
  if (params.get<INPAR::FLUID::LinearisationAction>("Linearisation")==INPAR::FLUID::Newton)
    is_newton_       = true;
  if (params.get<INPAR::FLUID::LinearisationAction>("Linearisation")==INPAR::FLUID::minimal)
    dserror("There is no LinearisationAction minimal in the fluid formulation");

  // set flags for formuation of the convective velocity term (conservative or convective)
  std::string convformstr = params.get<std::string>("form of convective term");
  if (convformstr =="conservative") is_conservative_ = true;

  // set flag for physical type of fluid flow
  physicaltype_ = params.get<INPAR::FLUID::PhysicalType>("Physical Type");
  if (((physicaltype_ != INPAR::FLUID::boussinesq) and (physicaltype_ != INPAR::FLUID::incompressible))
      and (is_stationary_ == true))
    dserror("physical type is not supported in stationary FLUID implementation.");

//----------------------------------------------------------------------
// get control parameters for time integration
//----------------------------------------------------------------------

  // get current time: n+alpha_F for generalized-alpha scheme, n+1 otherwise
  time_ = params.get<double>("total time",-1.0);

  // set global variable timefac to zero
  timefac_ = 0.0;

  if (is_stationary_ == false)
  {
    // get time-step length and time-integration parameters
    dt_                = params.get<double>("dt");
    const double theta = params.get<double>("theta",-1.0);
    omtheta_           = params.get<double>("omtheta",-1.0);

    // compute timefactor for left-hand side
    // in the case of BDF2 and generalized-alpha:
    // theta is set to the characteristic values in FLD::FluidImplicitTimeInt::PrepareTimeStep()

    // One-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = (alpha_F/alpha_M) * gamma * dt
    timefac_ = theta*dt_;

    if(is_genalpha_)
    {
      gamma_ =params.get<double>("gamma");
      alphaF_=params.get<double>("alphaF");
      alphaM_=params.get<double>("alphaM");
    }
    else
    {
      gamma_ =theta;
      alphaF_=1.0;
      alphaM_=1.0;
    }

    // in the case of not genalpha: afgdt = theta * dt_ = timefac_
    afgdt_=alphaF_*gamma_*dt_;
  }
  else
  {
    // timefac stationary = 1.0
    timefac_ = 1.0;
  }

  if (timefac_ < 0.0 or time_ < 0.0)
    dserror("Negative time-integration parameter or time-step length supplied");

// ---------------------------------------------------------------------
// get control parameters for stabilization and higher-order elements
//----------------------------------------------------------------------
  Teuchos::ParameterList& stablist = params.sublist("STABILIZATION");

  // no safety check necessary since all options are used
  tds_      = Teuchos::getIntegralValue<INPAR::FLUID::SubscalesTD>(stablist,"TDS");
  transient_= Teuchos::getIntegralValue<INPAR::FLUID::Transient>(stablist,"TRANSIENT");
  pspg_     = Teuchos::getIntegralValue<INPAR::FLUID::PSPG>(stablist,"PSPG");
  supg_     = Teuchos::getIntegralValue<INPAR::FLUID::SUPG>(stablist,"SUPG");
  vstab_    = Teuchos::getIntegralValue<INPAR::FLUID::VStab>(stablist,"VSTAB");
  cstab_    = Teuchos::getIntegralValue<INPAR::FLUID::CStab>(stablist,"CSTAB");
  cross_    = Teuchos::getIntegralValue<INPAR::FLUID::CrossStress>(stablist,"CROSS-STRESS");
  reynolds_ = Teuchos::getIntegralValue<INPAR::FLUID::ReynoldsStress>(stablist,"REYNOLDS-STRESS");

//-------------------------------
// get tau definition
//-------------------------------

  whichtau_ =  Teuchos::getIntegralValue<INPAR::FLUID::TauType>(stablist,"DEFINITION_TAU");
  // check if tau can be handled
  if (not(whichtau_ == INPAR::FLUID::tautype_franca_barrenechea_valentin_wall or
      INPAR::FLUID::tautype_franca_barrenechea_valentin_wall_wo_dt or
      INPAR::FLUID::tautype_bazilevs or
      INPAR::FLUID::tautype_bazilevs_wo_dt or
      INPAR::FLUID::tautype_franca_barrenechea_valentin_codina or
      INPAR::FLUID::tautype_oberai))
    dserror("Definition of Tau cannot be handled by the element");

  //stationary tau definition is switched on automatically
  if (is_stationary_ == true)
  {
    if (whichtau_ == INPAR::FLUID::tautype_franca_barrenechea_valentin_wall or
       whichtau_ == INPAR::FLUID::tautype_franca_barrenechea_valentin_wall_wo_dt)
    {
      whichtau_ = INPAR::FLUID::tautype_franca_barrenechea_valentin_wall_wo_dt;
    }
    else if (whichtau_ == INPAR::FLUID::tautype_bazilevs or
        whichtau_ == INPAR::FLUID::tautype_bazilevs_wo_dt)
    {
      whichtau_ = INPAR::FLUID::tautype_bazilevs_wo_dt;
    }
    else
      dserror("Actual tau definition (%i) is not compatible with the stationary option",
          whichtau_);
  }

  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (stablist.get<std::string>("STABTYPE") == "inconsistent") is_inconsistent_ = true;

  // set flags for potential evaluation of tau and material law at int. point
  // default value: evaluation at element center
  const std::string tauloc = stablist.get<std::string>("EVALUATION_TAU");
  if (tauloc == "integration_point") tau_gp_ = true;
  else                               tau_gp_ = false;
  const std::string matloc = stablist.get<std::string>("EVALUATION_MAT");
  if (matloc == "integration_point") mat_gp_ = true;
  else                               mat_gp_ = false;

//---------------------------------------------------------------------------------
// parameter for turbulence and subgrid-viscosity approach
//---------------------------------------------------------------------------------

  // get flag for fine-scale subgrid-viscosity approach
  {
    const std::string fssgvdef = params.get<std::string>("fs subgrid viscosity","No");

    if (fssgvdef == "Smagorinsky_all")        fssgv_ = INPAR::FLUID::smagorinsky_all;
    else if (fssgvdef == "Smagorinsky_small") fssgv_ = INPAR::FLUID::smagorinsky_small;
  }

  Teuchos::ParameterList& turbmodelparams    = params.sublist("TURBULENCE MODEL");

  if (fssgv_ != INPAR::FLUID::no_fssgv) Cs_ = turbmodelparams.get<double>("C_SMAGORINSKY",0.0);

  // the default action is no model
  turb_mod_action_ = INPAR::FLUID::no_model;

  // No turbulent flow: TURBULENCE_APPROACH = DNS
  if (turbmodelparams.get<std::string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    if (is_stationary_ == true)
      dserror("Stationary turbulent flow does not make any sense");

    // get Smagorinsky model parameter for fine-scale subgrid viscosity
    // (Since either all-scale Smagorinsky model (i.e., classical LES model
    // as will be inititalized below) or fine-scale Smagorinsky model is
    // used (and never both), the same input parameter can be exploited.)

    std::string& physical_turbulence_model = turbmodelparams.get<std::string>("PHYSICAL_MODEL");

    // --------------------------------------------------
    // standard constant coefficient Smagorinsky model
    if (physical_turbulence_model == "Smagorinsky")
    {
      // the classic Smagorinsky model only requires one constant parameter
      turb_mod_action_ = INPAR::FLUID::smagorinsky;
      Cs_              = turbmodelparams.get<double>("C_SMAGORINSKY");
    }
    // --------------------------------------------------
    // Smagorinsky model with van Driest damping
    else if (physical_turbulence_model == "Smagorinsky_with_van_Driest_damping")
    {
      // that's only implemented for turbulent channel flow
      // wall function length is hard implemented
      if (turbmodelparams.get<std::string>("CANONICAL_FLOW","no")
          !=
          "channel_flow_of_height_2")
          dserror("van_Driest_damping only for channel_flow_of_height_2\n");

      // for the Smagorinsky model with van Driest damping, we need
      // a viscous length to determine the y+ (heigth in wall units)
      turb_mod_action_ = INPAR::FLUID::smagorinsky_with_van_Driest_damping;

      // get parameters of model
      Cs_              = turbmodelparams.get<double>("C_SMAGORINSKY");
      l_tau_           = turbmodelparams.get<double>("CHANNEL_L_TAU");
    }

    // --------------------------------------------------
    // Smagorinsky model with dynamic Computation of Cs
    else if (physical_turbulence_model == "Dynamic_Smagorinsky")
    {
      turb_mod_action_ = INPAR::FLUID::dynamic_smagorinsky;

      // In the case of dynamic Smagorinsky:
      // Cs_ is calculated from Cs_sqrt_delta to compare it with the standard
      // it is stored in Cs_ after its calculation in CalcSubgrVisc
      Cs_ = 0.0;
    }
    else
    {
      dserror("Up to now, only Smagorinsky (constant coefficient with and without wall function as well as dynamic) is available");
    }
  } // end if(Classical LES)
}
