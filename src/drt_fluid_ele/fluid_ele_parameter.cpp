/*----------------------------------------------------------------------*/
/*!
\file fluid3_impl_parameter.cpp

\brief Evaluation of general fluid parameter

FluidEleParameter::SetParameter(Teuchos::ParameterList& params)
set all general fluid parameter once for all elements.

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#include <string>
#include <iostream>

#include "fluid_ele_parameter.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_io/io_pstream.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
std::map<int,Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> > DRT::ELEMENTS::FluidEleParameter::instances_;


//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> DRT::ELEMENTS::FluidEleParameter::Instance(int num)
{
  if (static_cast<int>(instances_.count(num))==0)
  {
    instances_.insert(std::pair<int,Teuchos::RCP<FluidEleParameter> >(num,Teuchos::rcp(new FluidEleParameter())));
  }
  return instances_.find(num)->second;
}



//----------------------------------------------------------------------*/
// private constructor of FluidEleParameter
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameter::FluidEleParameter()
  :
  set_general_fluid_parameter_(false),
  is_genalpha_(false),
  is_genalpha_np_(false),
  is_conservative_(false),
  is_stationary_(false),
  is_newton_(false),
  is_inconsistent_(false),
  reaction_(false),
  darcy_(false),
  reaction_topopt_(false),
  physicaltype_(INPAR::FLUID::incompressible),
  stabtype_(INPAR::FLUID::stabtype_nostab),
  tds_(INPAR::FLUID::subscales_quasistatic),
  transient_(INPAR::FLUID::inertia_stab_drop),
  pspg_(true),
  supg_(true),
  vstab_(INPAR::FLUID::viscous_stab_none),
  graddiv_(true),
  cross_(INPAR::FLUID::cross_stress_stab_none),
  reynolds_(INPAR::FLUID::reynolds_stress_stab_none),
  whichtau_(INPAR::FLUID::tau_not_defined),
  fssgv_(INPAR::FLUID::no_fssgv),
  viscreastabfac_(0.0),
  EOS_pres_(INPAR::FLUID::EOS_PRES_none),
  EOS_conv_stream_(INPAR::FLUID::EOS_CONV_STREAM_none),
  EOS_conv_cross_(INPAR::FLUID::EOS_CONV_CROSS_none),
  EOS_div_(INPAR::FLUID::EOS_DIV_none),
  EOS_whichtau_(INPAR::FLUID::EOS_tau_burman_fernandez),
  EOS_element_lenght_(INPAR::FLUID::EOS_he_max_dist_to_opp_surf),
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
  timefacpre_(1.0),
  turb_mod_action_(INPAR::FLUID::no_model),
  Cs_(0.0),
  Cs_averaged_(false),
  Ci_(0.0),
  include_Ci_(false),
  van_Driest_damping_(1.0),
  l_tau_(0.0),
  Cl_(0.0),
  Csgs_(0.0),
  Csgs_phi_(0.0),
  alpha_(0.0),
  CalcN_(false),
  N_(0.0),
  refvel_(INPAR::FLUID::strainrate),
  reflength_(INPAR::FLUID::cube_edge),
  c_nu_(1.0),
  c_diff_(1.0),
  near_wall_limit_(false),
  near_wall_limit_scatra_(false),
  B_gp_(false),
  beta_(0.0),
  mfs_is_conservative_(false),
  adapt_Csgs_phi_(false),
  meanCai_(0.0),
  update_mat_(false),
  conti_supg_(true),
  conti_cross_(INPAR::FLUID::cross_stress_stab_none),
  conti_reynolds_(INPAR::FLUID::reynolds_stress_stab_none),
  multifrac_loma_conti_(false),
  poro_conti_partint_(false)

{
}


//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementGeneralFluidParameter( Teuchos::ParameterList& params,
                                                                        int myrank )
{
  if(set_general_fluid_parameter_ == false)
    set_general_fluid_parameter_ = true;
  // For turbulent inflow generation,
  // this function is indeed two times called.
  // In this sepcial case, calling this function twice
  // is ok!
  else
  {
    if (myrank == 0)
      std::cout << std::endl << (" Warning: general fluid parameters should be set only once!!\n "
        "If you run a turbulent inflow generation, calling this function twice is ok!\n "
        "If you run a FPSI problem, calling this function twice is ok, too!\n"
        "Otherwise: Check, why you enter this function a second time!!!") << std::endl << std::endl;
  }

  //----------------------------------------------------------------------
  // get flags to switch on/off different fluid formulations
  //----------------------------------------------------------------------

  // set flag, time integration scheme
  timealgo_ = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "TimeIntegrationScheme");

  // set time integration scheme-specific element parameters
  if (timealgo_==INPAR::FLUID::timeint_stationary)
  {
    is_genalpha_ = false;
    is_stationary_ = true;
    is_genalpha_np_ = false;
  }
  else if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    is_genalpha_ = true;
    is_stationary_ = false;
    is_genalpha_np_ = false;
  }
  else if (timealgo_==INPAR::FLUID::timeint_npgenalpha)
  {
    is_genalpha_ = true;
    is_stationary_ = false;
    is_genalpha_np_ = true;
  }
  else
  {
    is_genalpha_ = false;
    is_stationary_ = false;
    is_genalpha_np_ = false;
  }

  // set flag for type of linearization (fixed-point-like or Newton)
  //std::string newtonstr   = params.get<std::string>("Linearisation");
  if (DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(params, "Linearisation")==INPAR::FLUID::Newton)
    is_newton_       = true;

  // set flags for formuation of the convective velocity term (conservative or convective)
  std::string convformstr = params.get<std::string>("form of convective term");
  if (convformstr =="conservative")
  {
    is_conservative_ = true;
    if(myrank==0)
    {
      std::cout << std::endl << "Warning: \n"
        "a) Using PSPG stabilization yields a conservative formulation (Hughes & Wells 2005)\n"
        "b) Instablities may occur for complex flow situations" << std::endl;
    }
  }

  // set flag for physical type of fluid flow
  physicaltype_ = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params, "Physical Type");
  if (
       (
            (physicaltype_ == INPAR::FLUID::loma)
         or (physicaltype_ == INPAR::FLUID::varying_density)
       )
      and (is_stationary_ == true)
     )
    dserror("physical type is not supported in stationary FLUID implementation.");

  if (is_genalpha_np_ and physicaltype_ == INPAR::FLUID::loma)
    dserror("the combination Np_Gen_Alpha and loma is not supported");

  if (not is_genalpha_ and physicaltype_ == INPAR::FLUID::loma)
    dserror("the combination OST and loma is said to be supported but does not work!!");

  if (is_genalpha_np_ and is_conservative_)
    dserror("the combination Np_Gen_Alpha and conservative flow is not supported");

  if (not is_stationary_ and is_conservative_ and physicaltype_ != INPAR::FLUID::incompressible)
  {
    if (myrank == 0)
     std::cout << std::endl << "Warning: missing time derivative terms in conservative formulation for variable density flows!" << std::endl;
  }

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------
  Teuchos::ParameterList stablist;
  if(params.get<int>("numfield")==0)
    stablist = params.sublist("RESIDUAL-BASED STABILIZATION");
  else if (params.get<int>("numfield")==1)
    stablist = params.sublist("POROUS-FLOW STABILIZATION");

  Teuchos::ParameterList& stablist_edgebased     = params.sublist("EDGE-BASED STABILIZATION");

  stabtype_ = DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(stablist, "STABTYPE");

  if (stabtype_ == INPAR::FLUID::stabtype_residualbased)
  {
    // no safety check necessary since all options are used
    tds_      = DRT::INPUT::IntegralValue<INPAR::FLUID::SubscalesTD>(stablist,"TDS");
    transient_= DRT::INPUT::IntegralValue<INPAR::FLUID::Transient>(stablist,"TRANSIENT");
    pspg_     = DRT::INPUT::IntegralValue<int>(stablist,"PSPG");
    supg_     = DRT::INPUT::IntegralValue<int>(stablist,"SUPG");
    vstab_    = DRT::INPUT::IntegralValue<INPAR::FLUID::VStab>(stablist,"VSTAB");
    rstab_    = DRT::INPUT::IntegralValue<INPAR::FLUID::RStab>(stablist,"RSTAB");
    graddiv_    = DRT::INPUT::IntegralValue<int>(stablist,"GRAD_DIV");
    cross_    = DRT::INPUT::IntegralValue<INPAR::FLUID::CrossStress>(stablist,"CROSS-STRESS");
    reynolds_ = DRT::INPUT::IntegralValue<INPAR::FLUID::ReynoldsStress>(stablist,"REYNOLDS-STRESS");

    // overrule higher_order_ele if input-parameter is set
    // this might be interesting for fast (but slightly
    // less accurate) computations
    is_inconsistent_ = DRT::INPUT::IntegralValue<int>(stablist,"INCONSISTENT");
    //-------------------------------
    // get tau definition
    //-------------------------------

    whichtau_ =  DRT::INPUT::IntegralValue<INPAR::FLUID::TauType>(stablist,"DEFINITION_TAU");
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
                         INPAR::FLUID::tau_codina_wo_dt or
                         INPAR::FLUID::tau_franca_madureira_valentin_badia_codina or
                         INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt))
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

    // get and check characteristic element length for stabilization parameter tau_Mu
    charelelengthu_ = DRT::INPUT::IntegralValue<INPAR::FLUID::CharEleLengthU>(stablist,"CHARELELENGTH_U");
    if (not(charelelengthu_ == INPAR::FLUID::streamlength_u or
                             INPAR::FLUID::volume_equivalent_diameter_u or
                             INPAR::FLUID::root_of_volume_u))
      dserror("Unknown characteristic element length for tau_Mu!");

    // get and check characteristic element length for stabilization parameter
    // tau_Mp and tau_C
    charelelengthpc_ = DRT::INPUT::IntegralValue<INPAR::FLUID::CharEleLengthPC>(stablist,"CHARELELENGTH_PC");
    if (not(charelelengthpc_ == INPAR::FLUID::streamlength_pc or
                             INPAR::FLUID::volume_equivalent_diameter_pc or
                             INPAR::FLUID::root_of_volume_pc))
      dserror("Unknown characteristic element length for tau_Mp and tau_C!");

    // in case of viscous and/or reactive stabilization, decide whether to use
    // GLS or USFEM and ensure compatibility of respective definitions
    if (vstab_ == INPAR::FLUID::viscous_stab_usfem or
        vstab_ == INPAR::FLUID::viscous_stab_usfem_only_rhs)
    {
      viscreastabfac_ = -1.0;
      if (rstab_ == INPAR::FLUID::reactive_stab_gls)
        dserror("inconsistent reactive and viscous stabilization!");
    }
    else if (vstab_ == INPAR::FLUID::viscous_stab_gls or
        vstab_ == INPAR::FLUID::viscous_stab_gls_only_rhs)
    {
      viscreastabfac_ = 1.0;
      if (rstab_ == INPAR::FLUID::reactive_stab_usfem)
        dserror("inconsistent reactive and viscous stabilization!");
    }
    else if (vstab_ == INPAR::FLUID::viscous_stab_none)
    {
      if (rstab_ == INPAR::FLUID::reactive_stab_usfem)    viscreastabfac_ = -1.0;
      else if (rstab_ == INPAR::FLUID::reactive_stab_gls) viscreastabfac_ =  1.0;
    }

    // case of xfem check whether additional xfem-stabilization terms in the form of
    // edge-based terms are activated (i.e., ghost penalties)
    if (stablist_edgebased.get<std::string>("EOS_PRES") == "xfem_gp" or
        stablist_edgebased.get<std::string>("EOS_CONV_STREAM") == "xfem_gp" or
        stablist_edgebased.get<std::string>("EOS_CONV_CROSS") == "xfem_gp" or
        stablist_edgebased.get<std::string>("EOS_DIV") == "xfem_gp")
    {
        EOS_pres_         = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Pres>(stablist_edgebased,"EOS_PRES");
        EOS_conv_stream_  = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Conv_Stream>(stablist_edgebased,"EOS_CONV_STREAM");
        EOS_conv_cross_   = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Conv_Cross>(stablist_edgebased,"EOS_CONV_CROSS");
        EOS_div_          = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Div>(stablist_edgebased,"EOS_DIV");

        EOS_whichtau_       = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_TauType>(stablist_edgebased, "EOS_DEFINITION_TAU");
    }
    // setting the EOS element length outside of the if-statement is not very beautiful,
    // but there is an input parameter in the XFEM STABILIZATION section, which requires that
    // this parameter has been set
    EOS_element_lenght_ = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_ElementLength>(stablist_edgebased, "EOS_H_DEFINITION");
  }
  else if (stabtype_ == INPAR::FLUID::stabtype_edgebased)
  {
    if (myrank==0)
    {
      IO::cout << "+----------------------------------------------------------------------------------+\n";
      IO::cout << " Edge-based stabilization: all residual-based stabilization terms are switched off!\n";
      IO::cout << "+----------------------------------------------------------------------------------+\n" << IO::endl;
    }
    //---------------------------------
    // if edge-based stabilization is selected, all residual-based stabilization terms
    // are switched off
    pspg_ = false;
    supg_ = false;
    vstab_ = INPAR::FLUID::viscous_stab_none;
    rstab_ = INPAR::FLUID::reactive_stab_none;
    graddiv_ = false;
    cross_ = INPAR::FLUID::cross_stress_stab_none;
    reynolds_ = INPAR::FLUID::reynolds_stress_stab_none;
    tds_ = INPAR::FLUID::subscales_quasistatic;
    transient_ = INPAR::FLUID::inertia_stab_drop;
    is_inconsistent_ = false;

    // --------------------------------
    // edge-based fluid stabilization can be used as standard fluid stabilization or
    // as ghost-penalty stabilization in addition to residual-based stabilizations in the XFEM

    // set parameters if single stabilization terms are switched on/off or which type of stabilization is chosen
    EOS_pres_         = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Pres>(stablist_edgebased,"EOS_PRES");
    EOS_conv_stream_  = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Conv_Stream>(stablist_edgebased,"EOS_CONV_STREAM");
    EOS_conv_cross_   = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Conv_Cross>(stablist_edgebased,"EOS_CONV_CROSS");
    EOS_div_          = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Div>(stablist_edgebased,"EOS_DIV");

    EOS_element_lenght_ = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_ElementLength>(stablist_edgebased, "EOS_H_DEFINITION");
    EOS_whichtau_       = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_TauType>(stablist_edgebased, "EOS_DEFINITION_TAU");

  }
  else if (stabtype_ == INPAR::FLUID::stabtype_nostab)
  {
    if (myrank==0)
    {
      IO::cout << "+----------------------------------------------------------------------------------+\n";
      IO::cout << "+                                   WARNING\n";
      IO::cout << " No stabilization selected: all stabilization terms are switched off!\n";
      IO::cout << "                            BACI says: Good luck!\n";
      IO::cout << "+----------------------------------------------------------------------------------+\n" << IO::endl;
    }
    pspg_ = false;
    supg_ = false;
    vstab_ = INPAR::FLUID::viscous_stab_none;
    rstab_ = INPAR::FLUID::reactive_stab_none;
    graddiv_ = false;
    cross_ = INPAR::FLUID::cross_stress_stab_none;
    reynolds_ = INPAR::FLUID::reynolds_stress_stab_none;
    tds_ = INPAR::FLUID::subscales_quasistatic;
    transient_ = INPAR::FLUID::inertia_stab_drop;
    is_inconsistent_ = false;
  }
  else
   dserror("Unknown stabilization type");


  //---------------------------------
  // safety checks for time-dependent subgrid scales
  if ((tds_ == INPAR::FLUID::subscales_time_dependent) or (transient_ != INPAR::FLUID::inertia_stab_drop))
  {
    if (not is_genalpha_np_)
      dserror("time dependent subscales does not work for OST/AfGenAlpha/BDF2/Stationary. \nOne need to look for bugs");
  }


  //---------------------------------
  // set flags for potential evaluation of tau and material law at int. point
  // default value: evaluation at element center
  const std::string tauloc = stablist.get<std::string>("EVALUATION_TAU");
  if (tauloc == "integration_point") tau_gp_ = true;
  else                               tau_gp_ = false;
  const std::string matloc = stablist.get<std::string>("EVALUATION_MAT");
  if (matloc == "integration_point") mat_gp_ = true;
  else                               mat_gp_ = false;


}

void DRT::ELEMENTS::FluidEleParameter::SetElementTimeParameter( Teuchos::ParameterList& params )
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

    //-----------------------------------------------------------------------
    //       |          timefac         |  timefacpre     |    timefacrhs   |
    // ----------------------------------------------------------------------
    // OST   |                        dt*theta                              |
    //-----------------------------------------------------------------------
    // BDF2  |                        2/3 * dt                              |
    //-----------------------------------------------------------------------
    // Af GA |          alphaF*gamma*dt/alphaM            | gamma*dt/alphaM |
    //----------------------------------------------------------------------
    // NP GA | alphaF*gamma*dt/alphaM   | gamma*dt/alphaM | gamma*dt/alphaM |
    //-----------------------------------------------------------------------

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

    // timeint_gen_alpha = p(n+1) (Peter's genalpha)
    if (timealgo_ == INPAR::FLUID::timeint_npgenalpha)
    {
      // if not generalized-alpha: timefacrhs_=theta * dt_ = timefac_
      timefacpre_ = gamma_/alphaM_*dt_;
      timefacrhs_ = gamma_/alphaM_*dt_;
    }
    else if(timealgo_ == INPAR::FLUID::timeint_afgenalpha)
    {
      timefacpre_ = gamma_*alphaF_/alphaM_*dt_;
      timefacrhs_ = gamma_/alphaM_*dt_;
    }
    else
    {
      // if not generalized-alpha: timefacmat_p_=theta * dt_ = timefac_
      timefacpre_ = gamma_*alphaF_/alphaM_*dt_;
      // if not generalized-alpha: timefacrhs_=theta * dt_ = timefac_
      timefacrhs_ = gamma_*alphaF_/alphaM_*dt_;
    }
  }
  else // is_stationary == true
  {
    // set timefactor for stationary case to 1.0
    timefac_ = 1.0;
    timefacrhs_ = 1.0;
  }

  if (dt_ < 0.0 or theta_ < 0.0 or time_ < 0.0 or omtheta_ < 0.0 or gamma_ < 0.0
      or alphaF_ < 0.0 or alphaM_ < 0.0)
  {
    std::cout<<"dt_: "<<dt_<<std::endl;
    std::cout<<"theta_ "<<theta_<<std::endl;
    std::cout<<"time_ "<<time_<<std::endl;
    std::cout<<"omtheta_ "<<omtheta_<<std::endl;
    std::cout<<"gamma_ "<<gamma_<<std::endl;
    std::cout<<"alphaF_ "<<alphaF_<<std::endl;
    std::cout<<"alphaM_ "<<alphaM_<<std::endl;
    dserror("Negative (or no) time-integration parameter or time-step length supplied");
  }
}

//----------------------------------------------------------------------*
//  set turbulence parameters                            rasthofer 11/11|
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementTurbulenceParameter( Teuchos::ParameterList& params )
{
  // get parameter lists
  Teuchos::ParameterList& turbmodelparams = params.sublist("TURBULENCE MODEL");
  Teuchos::ParameterList& turbmodelparamssgvisc = params.sublist("SUBGRID VISCOSITY");
  Teuchos::ParameterList& turbmodelparamsmfs = params.sublist("MULTIFRACTAL SUBGRID SCALES");

  //---------------------------------------------------------------------------------
  // parameter for subgrid-viscosity approach
  //---------------------------------------------------------------------------------

  // get flag for fine-scale subgrid-viscosity approach
  {
    const std::string fssgvdef = turbmodelparams.get<std::string>("FSSUGRVISC","No");

    if (fssgvdef == "Smagorinsky_all")        fssgv_ = INPAR::FLUID::smagorinsky_all;
    else if (fssgvdef == "Smagorinsky_small") fssgv_ = INPAR::FLUID::smagorinsky_small;
  }

  // get Smagorinsky model parameter for fine-scale subgrid viscosity
  // (Since either all-scale Smagorinsky model (i.e., classical LES model
  // as will be inititalized below) or fine-scale Smagorinsky model is
  // used (and never both), the same input parameter can be exploited.)
  if (fssgv_ != INPAR::FLUID::no_fssgv) Cs_ = turbmodelparamssgvisc.get<double>("C_SMAGORINSKY",0.0);

  //---------------------------------------------------------------------------------
  // parameter for turbulence approach
  //---------------------------------------------------------------------------------

  // the default action is no model
  turb_mod_action_ = INPAR::FLUID::no_model;

  // No turbulent flow: TURBULENCE_APPROACH = DNS
  if (turbmodelparams.get<std::string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    if (is_stationary_ == true)
      dserror("Stationary turbulent flow does not make any sense");

    std::string& physical_turbulence_model = turbmodelparams.get<std::string>("PHYSICAL_MODEL");

    // --------------------------------------------------
    // standard constant coefficient Smagorinsky model
    if (physical_turbulence_model == "Smagorinsky")
    {
      // the classic Smagorinsky model only requires one constant parameter
      turb_mod_action_ = INPAR::FLUID::smagorinsky;
      Cs_              = turbmodelparamssgvisc.get<double>("C_SMAGORINSKY");
      include_Ci_ = DRT::INPUT::IntegralValue<int>(turbmodelparamssgvisc,"C_INCLUDE_CI");
      Ci_              = turbmodelparamssgvisc.get<double>("C_YOSHIZAWA");
    }
    // --------------------------------------------------
    // Smagorinsky model with van Driest damping
    else if (physical_turbulence_model == "Smagorinsky_with_van_Driest_damping")
    {
      // that's only implemented for turbulent channel flow
      // wall function length is hard implemented
      if (turbmodelparamssgvisc.get<std::string>("CANONICAL_FLOW","no")
          !=
          "channel_flow_of_height_2")
          dserror("van_Driest_damping only for channel_flow_of_height_2\n");

      // for the Smagorinsky model with van Driest damping, we need
      // a viscous length to determine the y+ (heigth in wall units)
      turb_mod_action_ = INPAR::FLUID::smagorinsky_with_van_Driest_damping;

      // get parameters of model
      Cs_              = turbmodelparamssgvisc.get<double>("C_SMAGORINSKY");
      l_tau_           = turbmodelparamssgvisc.get<double>("CHANNEL_L_TAU");
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
      Cs_averaged_ = DRT::INPUT::IntegralValue<int>(turbmodelparamssgvisc,"C_SMAGORINSKY_AVERAGED");
      Ci_ = turbmodelparamssgvisc.get<double>("C_YOSHIZAWA");
      include_Ci_ = DRT::INPUT::IntegralValue<int>(turbmodelparamssgvisc,"C_INCLUDE_CI");
    }
    else if (physical_turbulence_model == "Scale_Similarity")
    {
      turb_mod_action_ = INPAR::FLUID::scale_similarity;
      Cl_ = turbmodelparamsmfs.get<double>("C_SCALE_SIMILARITY");
    }
    else if (physical_turbulence_model == "Scale_Similarity_basic")
    {
      turb_mod_action_ = INPAR::FLUID::scale_similarity_basic;
      Cl_ = turbmodelparamsmfs.get<double>("C_SCALE_SIMILARITY");
    }
    else if (physical_turbulence_model == "Multifractal_Subgrid_Scales")
    {
      turb_mod_action_ = INPAR::FLUID::multifractal_subgrid_scales;

      // get parameters of model
      Csgs_ = turbmodelparamsmfs.get<double>("CSGS");
      Csgs_phi_ = turbmodelparamsmfs.get<double>("CSGS_PHI");
      adapt_Csgs_phi_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs,"ADAPT_CSGS_PHI");

      if (turbmodelparamsmfs.get<std::string>("SCALE_SEPARATION") == "algebraic_multigrid_operator")
       alpha_ = 3.0;
      else if (turbmodelparamsmfs.get<std::string>("SCALE_SEPARATION") == "box_filter"
            or turbmodelparamsmfs.get<std::string>("SCALE_SEPARATION") == "geometric_multigrid_operator")
       alpha_ = 2.0;
      else
       dserror("Unknown filter type!");

      CalcN_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs,"CALC_N");

      N_ = turbmodelparamsmfs.get<double>("N");

      if (turbmodelparamsmfs.get<std::string>("REF_VELOCITY") == "strainrate")
       refvel_ = INPAR::FLUID::strainrate;
      else if (turbmodelparamsmfs.get<std::string>("REF_VELOCITY") == "resolved")
       refvel_ = INPAR::FLUID::resolved;
      else if (turbmodelparamsmfs.get<std::string>("REF_VELOCITY") == "fine_scale")
       refvel_ = INPAR::FLUID::fine_scale;
      else
       dserror("Unknown velocity!");

      if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "cube_edge")
       reflength_ = INPAR::FLUID::cube_edge;
      else if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "sphere_diameter")
       reflength_ = INPAR::FLUID::sphere_diameter;
      else if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "streamlength")
       reflength_ = INPAR::FLUID::streamlength;
      else if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "gradient_based")
       reflength_ = INPAR::FLUID::gradient_based;
      else if (turbmodelparamsmfs.get<std::string>("REF_LENGTH") == "metric_tensor")
       reflength_ = INPAR::FLUID::metric_tensor;
      else
       dserror("Unknown length!");

      c_nu_ = turbmodelparamsmfs.get<double>("C_NU");
      c_diff_ = turbmodelparamsmfs.get<double>("C_DIFF"); //loma only

      near_wall_limit_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs,"NEAR_WALL_LIMIT");
      near_wall_limit_scatra_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs,"NEAR_WALL_LIMIT_CSGS_PHI");

      if (turbmodelparamsmfs.get<std::string>("EVALUATION_B") == "element_center")
      B_gp_ = false;
      else if (turbmodelparamsmfs.get<std::string>("EVALUATION_B") == "integration_point")
      B_gp_ = true;
      else
        dserror("Unknown evaluation point!");

      beta_ = turbmodelparamsmfs.get<double>("BETA");

      if (turbmodelparamsmfs.get<std::string>("CONVFORM") == "conservative")
       mfs_is_conservative_ = true;
      else
       mfs_is_conservative_ = false;
    }
    else
    {
      dserror("Up to now, only Smagorinsky, Scale Similarity and Multifractal Subgrid Scales are available");
    }
  } // end if(Classical LES)
}


//----------------------------------------------------------------------*
//  set loma parameters                                  rasthofer 03/12|
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementLomaParameter( Teuchos::ParameterList& params )
{
  // get parameter lists
  Teuchos::ParameterList& lomaparams = params.sublist("LOMA");
  Teuchos::ParameterList& stabparams = params.sublist("RESIDUAL-BASED STABILIZATION");
  Teuchos::ParameterList& turbmodelparamsmfs = params.sublist("MULTIFRACTAL SUBGRID SCALES");

  //---------------------------------------------------------------------------------
  // material update with subgrid-scale temperature
  //---------------------------------------------------------------------------------

  update_mat_ = lomaparams.get<bool>("update material",false);

  //---------------------------------------------------------------------------------
  // parameter for additional rbvmm terms in continuity equation
  //---------------------------------------------------------------------------------

  conti_supg_     = DRT::INPUT::IntegralValue<int>(stabparams,"LOMA_CONTI_SUPG");
  conti_cross_    = DRT::INPUT::IntegralValue<INPAR::FLUID::CrossStress>(stabparams,"LOMA_CONTI_CROSS_STRESS");
  conti_reynolds_ = DRT::INPUT::IntegralValue<INPAR::FLUID::ReynoldsStress>(stabparams,"LOMA_CONTI_REYNOLDS_STRESS");

  //---------------------------------------------------------------------------------
  // parameter for additional multifractal subgrid-scale terms
  //---------------------------------------------------------------------------------

  if (turb_mod_action_ == INPAR::FLUID::multifractal_subgrid_scales)
   multifrac_loma_conti_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs,"LOMA_CONTI");

  return;
}


//----------------------------------------------------------------------*
//  set poro parameters                                  vuong 11/12|
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementPoroParameter( Teuchos::ParameterList& params )
{
  poro_conti_partint_ = params.get<bool>("conti partial integration",false);
  reaction_= true;
  reaction_topopt_= false;
  darcy_= true;
  graddiv_=false;

  return;
}


//----------------------------------------------------------------------*
//  set topopt parameters                               winklmaier 07/13|
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementTopoptParameter( Teuchos::ParameterList& params )
{
  topopt_params_[0] = params.get<double>("MIN_PORO");
  topopt_params_[1] = params.get<double>("MAX_PORO");
  topopt_params_[2] = params.get<double>("SMEAR_FAC");
  reaction_= true;
  reaction_topopt_= true;
  darcy_= false;

  return;
}


//----------------------------------------------------------------------*/
// print fluid parameter to screen (AE 01-11)
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::PrintFluidParameter()
{
    std::cout << std::endl << "|-----------------------------------------------------------------------------" << std::endl;
    std::cout << "|  General Fluid parameter: " << std::endl;
    std::cout << "|-----------------------------------------------------------------------------" << std::endl;
    //! flag SetGeneralParameter was called
    std::cout << "|    method SetElmentGeneralFluidParameter was called:    " << set_general_fluid_parameter_ << std::endl;
    //! flag to (de)activate generalized-alpha time-integration scheme
    std::cout << "|    generalized alpha time integration active:    " << is_genalpha_ << std::endl;
    //! flag to (de)activate generalized-alpha-np1 time-integration scheme
    std::cout << "|    generalized alpha time integration np:    " << is_genalpha_np_ << std::endl;
    //! flag to (de)activate conservative formulation
    std::cout << "|    conservative formulation:    " << is_conservative_ << std::endl;
    //! flag to (de)activate stationary formulation
    std::cout << "|    steady state:    " << is_stationary_ << std::endl;
    //! flag to (de)activate Newton linearization
    std::cout << "|    Newton linearization:    " << is_newton_ << std::endl;
    //! flag to (de)activate second derivatives
    std::cout << "|    use inconsistent:    " << is_inconsistent_ << std::endl;
    //! flag to (de)activate potential reactive terms
    std::cout << "|    reaction term:    " << reaction_ << std::endl;
    //! flag to (de)aktivate porous darcy flow
    std::cout << "|    darcy equation:    " << darcy_ << std::endl;
    //! flag to (de)aktivate reaction due to topology optimization
    std::cout << "|    reaction term due to topology optimization:    " << reaction_topopt_ << std::endl;
    //! matrix with values for computation of porosity with respect to topopt density
    std::cout << "|    topology optimization porosity parameter: " << topopt_params_[0]
        << " " << topopt_params_[1] << " " << topopt_params_[2] << std::endl;
    //! Flag for physical type of the fluid flow (incompressible, loma, varying_density, Boussinesq)
    std::cout << "|    physical type:    "<< physicaltype_ << std::endl;
    //! Flag to (de)activate time-dependent subgrid stabilization
    std::cout << "|    time-dependent subgrid stabilization:    " << tds_ << std::endl;
    //! Flag to (de)activate time-dependent term in large-scale momentum equation
    std::cout << "|    time dependent term:    " << transient_ << std::endl;
    //! Flag to (de)activate PSPG stabilization
    std::cout << "|    PSPG:    " << pspg_ << std::endl;
    //! Flag to (de)activate SUPG stabilization
    std::cout << "|    SUPG:    " << supg_<< std::endl ;
    //! Flag to (de)activate viscous term in residual-based stabilization
    std::cout << "|    VSTAB:    " << vstab_ << std::endl;
    //! Flag to (de)activate least-squares stabilization of continuity equation
    std::cout << "|    Grad-Div-Stab:    " << graddiv_ << std::endl ;
    //! Flag to (de)activate reactive term in residual-based stabilization
    std::cout << "|    reactive stabilization:    " << rstab_ << std::endl;
    //! Flag to (de)activate cross-stress term -> residual-based VMM
    std::cout << "|    cross-stress term:    " << cross_ << std::endl;
    //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
    std::cout << "|    Reynolds-stress term:    " << reynolds_ << std::endl;
    //! (sign) factor for viscous and reactive stabilization terms
    std::cout << "|    viscous and reactive stabilization factor:    " << viscreastabfac_ << std::endl;
    //! Flag to define tau
    std::cout << "|    Definition of stabilization parameter:    " << whichtau_ << std::endl;
    //! flag to (de)activate fine-scale subgrid viscosity
    std::cout << "|    fine-scale subgrid viscosity::    " << fssgv_ << std::endl;
    //! flag for material evaluation at Gaussian integration points
    std::cout << "|    material evaluation at Gaussian integration points:    " << mat_gp_ << std::endl;
    //! flag for stabilization parameter evaluation at Gaussian integration points
    std::cout << "|    stabilization parameter evaluation at Gaussian integration points:    " << tau_gp_ << std::endl;
    std::cout << "|---------------------------------------------------------------------------" << std::endl;

    std::cout << std::endl << "|---------------------------------------------------------------------------" << std::endl;
    std::cout << "|  Time parameter: " << std::endl;
    std::cout << "|---------------------------------------------------------------------------" << std::endl;
    //! time algorithm
    std::cout << "|    time algorithm:    " << timealgo_ << std::endl;
    //! actual time to evaluate the body BC
    std::cout << "|    time:    " << time_ << std::endl;
    //! time-step length
    std::cout << "|    time step:    " << dt_ << std::endl;
    //! timefac = dt_ * ("pseudo"-)theta_
    std::cout << "|    time factor:    " << timefac_ << std::endl;
    //! factor for left-hand side due to one-step-theta time-integration scheme
    std::cout << "|    theta:    " << theta_ << std::endl;
    //! factor for right-hand side due to one-step-theta time-integration scheme
    std::cout << "|    (1-theta):    " << omtheta_ << std::endl;
    //! generalised-alpha parameter (connecting velocity and acceleration)
    std::cout << "|    gamma:    " << gamma_ << std::endl;
    //! generalised-alpha parameter (velocity)
    std::cout << "|    alpha_F:    " << alphaF_ << std::endl;
    //! generalised-alpha parameter (acceleration)
    std::cout << "|    alpha_M:    " << alphaM_ << std::endl;
    //! generalised-alpha parameter, alphaF_*gamma_*dt_
    std::cout << "|    time factor mat_u:    " << afgdt_ << std::endl;
    //! time integration factor for the right hand side (boundary elements)
    std::cout << "|    time factor rhs:    " << timefacrhs_ << std::endl;
    //! time integration factor for the left hand side (pressure)
    std::cout << "|    time factor mat_p:    " << timefacpre_ << std::endl;
    std::cout << "|---------------------------------------------------------------------------" << std::endl;

    std::cout << std::endl << "|---------------------------------------------------------------------------" << std::endl;
    std::cout << "|  Turbulence parameter: " << std::endl;
    std::cout << "|---------------------------------------------------------------------------" << std::endl;
    //! flag to define turbulence model
    std::cout << "|    turbulence model:    " << turb_mod_action_ << std::endl;
    //! smagorinsky constant
    std::cout << "|    smagorinsky constant:    " << Cs_ << std::endl;
    //! comment missing
    std::cout << "|    Cs_averaged_ is    " << Cs_averaged_ << std::endl;
    //! scale similarity constant
    std::cout << "|    Cl_ is    " << Cl_ << std::endl;
    /// multifractal subgrid-scales
    std::cout << "|    Csgs_ is    " << Csgs_ << std::endl;
    //! comment missing
    std::cout << "|    alpha_ is    " << alpha_ << std::endl;
    //! comment missing
    std::cout << "|    CalcN_ is    " << CalcN_ << std::endl;
    //! comment missing
    std::cout << "|    N_ is    " << N_ << std::endl;
    //! comment missing
    std::cout << "|    refvel_ is    " << refvel_ << std::endl;
    //! comment missing
    std::cout << "|    reflength_ is    " << reflength_ << std::endl;
    //! comment missing
    std::cout << "|    c_nu_ is    " << c_nu_ << std::endl;
    //! comment missing
    std::cout << "|    near_wall_limit_ is    " << near_wall_limit_ << std::endl;
    //! comment missing
    std::cout << "|    B_gp_ is    " << B_gp_ << std::endl;
    //! comment missing
    std::cout << "|    beta_ is    " << beta_ << std::endl;
    //! channel length to normalize the normal wall distance
    std::cout << "|    channel length to normalize the normal wall distance:    " << l_tau_ << std::endl;
    std::cout << "|---------------------------------------------------------------------------" << std::endl;

}
