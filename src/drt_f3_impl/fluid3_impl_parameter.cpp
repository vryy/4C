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
#include <iostream>

#include "fluid3_impl_parameter.H"
#include "../drt_lib/drt_dserror.H"

using namespace std;
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
  set_general_fluid_parameter_(false),
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
  time_(-1.0),
  dt_(0.0),
  timefac_(0.0),
  theta_(0.0),
  omtheta_(0.0),
  gamma_(0.0),
  alphaF_(0.0),
  alphaM_(0.0),
  afgdt_(1.0),
  timefacrhs_(1.0),
  timefacmat_p_(1.0),
  Cs_(0.0),
  l_tau_(0.0)
{
}


//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3ImplParameter::SetElementGeneralFluidParameter( Teuchos::ParameterList& params )
{
  // routine SetGeneralFluidParameter was called once
  if(set_general_fluid_parameter_ == false)
    set_general_fluid_parameter_ = true;
  else
    dserror("SetGeneralFluidParameter() is supposed to be called only once during runtime");

//----------------------------------------------------------------------
// get flags to switch on/off different fluid formulations
//----------------------------------------------------------------------

  // set flag, time integration scheme
  timealgo_ = params.get<INPAR::FLUID::TimeIntegrationScheme>("TimeIntegrationScheme");

  // set time integration scheme-specific element parameters
  if (timealgo_==INPAR::FLUID::timeint_stationary)
  {
    is_genalpha_ = false;
    is_stationary_ = true;
  }
  else if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    is_genalpha_ = true;
    is_stationary_ = false;
  }
  else if (timealgo_==INPAR::FLUID::timeint_gen_alpha)
  {
    is_genalpha_ = true;
    is_stationary_ = false;
  }
  else
  {
    is_genalpha_ = false;
    is_stationary_ = false;
  }

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
  if (not(whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins or
                       INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt or
                       INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen or
                       INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt or
                       INPAR::FLUID::tau_taylor_hughes_zarins_scaled or
                       INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt or
                       INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall or
                       INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt or
                       INPAR::FLUID::tau_shakib_hughes_codina or
                       INPAR::FLUID::tau_shakib_hughes_codina_wo_dt or
                       INPAR::FLUID::tau_codina or
                       INPAR::FLUID::tau_codina_wo_dt))
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
    else if (physical_turbulence_model == "Scale_Similarity")
    {
      turb_mod_action_ = INPAR::FLUID::scale_similarity;
      Cl_ = turbmodelparams.get<double>("C_SCALE_SIMILARITY");
    }
    else if (physical_turbulence_model == "Mixed_Scale_Similarity_Eddy_Viscosity_Model")
    {
      turb_mod_action_ = INPAR::FLUID::mixed_scale_similarity_eddy_viscosity_model;
      Cl_ = turbmodelparams.get<double>("C_SCALE_SIMILARITY");
      Cs_ = turbmodelparams.get<double>("C_SMAGORINSKY");
    }
    else
    {
      dserror("Up to now, only Smagorinsky (constant coefficient with and without wall function as well as dynamic) and Scale Similarity (also combined with Smagorinsky) are available");
    }
  } // end if(Classical LES)
}

void DRT::ELEMENTS::Fluid3ImplParameter::SetElementTimeParameter( Teuchos::ParameterList& params )
{
  // second check: timealgo
  // work around to use SetTimeParameter in GenaAlpha (Neumann BC)
  if(set_general_fluid_parameter_!= true)
    dserror("General fluid parameter are not set yet!!");

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------

  // get current time: n+alpha_F for generalized-alpha scheme, n+1 otherwise
  time_ = params.get<double>("total time",-1.0);

  // set global variable timefac to zero
  timefac_ = 0.0;

  if (not is_stationary_)
  {
    // get time-step length and time-integration parameters
    dt_      = params.get<double>("dt",-1.0);
    theta_   = params.get<double>("theta",-1.0);
    omtheta_ = params.get<double>("omtheta",-1.0);

    // compute timefactor for left-hand side:
    // one-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = (alpha_F/alpha_M) * gamma * dt
    // (For BDF2 and generalized-alpha, theta was already computed
    //  accordingly in FLD::FluidImplicitTimeInt::PrepareTimeStep().)
    timefac_ = theta_*dt_;

    // compute generalized-alpha-related values and set them appropriately
    // otherwise
    if (is_genalpha_)
    {
      gamma_  = params.get<double>("gamma",-1.0);
      alphaF_ = params.get<double>("alphaF",-1.0);
      alphaM_ = params.get<double>("alphaM",-1.0);
    }
    else
    {
      gamma_  = theta_;
      alphaF_ = 1.0;
      alphaM_ = 1.0;
    }

    // if not generalized-alpha: afgdt = theta * dt_ = timefac_
    // Peter's generalized alpha: timefacmat_u_ for velocity terms
    afgdt_=alphaF_*gamma_*dt_;

    if (timealgo_ == INPAR::FLUID::timeint_gen_alpha)
    {
      // if not generalized-alpha: timefacrhs_=theta * dt_ = timefac_
      timefacrhs_ = 1.0;
      timefacmat_p_ = gamma_*dt_;
    }
    else
    {
      // if not generalized-alpha: timefacrhs_=theta * dt_ = timefac_
      timefacrhs_ = gamma_/alphaM_*dt_;
      // if not generalized-alpha: timefacmat_p_=theta * dt_ = timefac_
      timefacmat_p_ = alphaF_*gamma_/alphaM_*dt_;
    }
  }
  else // is_stationary == true
  {
    // set timefactor for stationary case to 1.0
    timefac_ = 1.0;
  }

  if (dt_ < 0.0 or theta_ < 0.0 or time_ < 0.0 or omtheta_ < 0.0 or gamma_ < 0.0
      or alphaF_ < 0.0 or alphaM_ < 0.0)
    dserror("Negative (or no) time-integration parameter or time-step length supplied");
}
