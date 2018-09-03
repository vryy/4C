/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter.cpp

\brief Setting of general fluid parameter for element evaluation

This file has to contain all parameters called in fluid_ele_calc.cpp.
Additional parameters required in derived classes of FluidEleCalc have to
be set in problem specific parameter lists derived from this class.

\level 2

<pre>

\maintainer Volker Gravemeier
             vgravem@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#include <string>
#include <iostream>
#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io_pstream.H"
#include "fluid_ele_parameter.H"


//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameter::FluidEleParameter()
    : set_general_fluid_parameter_(false),
      physicaltype_(INPAR::FLUID::incompressible),
      stabtype_(INPAR::FLUID::stabtype_nostab),  // stabilization parameters
      is_conservative_(false),
      is_newton_(false),
      is_inconsistent_(false),
      reaction_(false),
      oseenfieldfuncno_(-1),
      is_reconstructder_(false),
      tds_(INPAR::FLUID::subscales_quasistatic),
      transient_(INPAR::FLUID::inertia_stab_drop),
      pspg_(true),
      supg_(true),
      vstab_(INPAR::FLUID::viscous_stab_none),
      rstab_(INPAR::FLUID::reactive_stab_none),
      graddiv_(true),
      cross_(INPAR::FLUID::cross_stress_stab_none),
      reynolds_(INPAR::FLUID::reynolds_stress_stab_none),
      whichtau_(INPAR::FLUID::tau_not_defined),
      charelelengthu_(INPAR::FLUID::streamlength_u),
      charelelengthpc_(INPAR::FLUID::volume_equivalent_diameter_pc),
      viscreastabfac_(0.0),
      ppp_(false),
      mat_gp_(false),             // standard evaluation of the material at the element center
      tau_gp_(false),             // standard evaluation of tau at the element center
      interface_thickness_(0.0),  // two phase parameters
      enhanced_gaussrule_(false),
      include_surface_tension_(false),  // include the surface tension in the calculations.
      eval_surfacetension_(INPAR::TWOPHASE::surface_tension_approx_none),
      turb_mod_action_(INPAR::FLUID::no_model),  // turbulence parameters
      Cs_(0.0),
      Cs_averaged_(false),
      Ci_(0.0),
      include_Ci_(false),
      van_Driest_damping_(1.0),
      l_tau_(0.0),
      fssgv_(INPAR::FLUID::no_fssgv),
      vrfi_(INPAR::FLUID::cuberootvol),
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
      consistent_mfs_residual_(false),
      update_mat_(false),
      conti_supg_(true),
      conti_cross_(INPAR::FLUID::cross_stress_stab_none),
      conti_reynolds_(INPAR::FLUID::reynolds_stress_stab_none),
      multifrac_loma_conti_(false),
      reaction_topopt_(false)
{
  // we have to know the time parameters here to check for illegal combinations
  fldparatimint_ = DRT::ELEMENTS::FluidEleParameterTimInt::Instance();
}

//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementGeneralFluidParameter(
    Teuchos::ParameterList& params, int myrank)
{
  if (set_general_fluid_parameter_ == false) set_general_fluid_parameter_ = true;
  // For turbulent inflow generation,
  // this function is indeed two times called.
  // In this sepcial case, calling this function twice
  // is ok!
  else
  {
    if (myrank == 0)
      std::cout
          << std::endl
          << (" Warning: general fluid parameters should be set only once!!\n "
              " If you run a turbulent inflow generation, calling this function twice is ok!\n ")
          << std::endl
          << std::endl;
  }

  // set flags for formulation of the convective velocity term (conservative or convective)
  std::string convformstr = params.get<std::string>("form of convective term");
  if (convformstr == "conservative")
  {
    is_conservative_ = true;
    if (myrank == 0)
    {
      std::cout
          << std::endl
          << "Warning: \n"
             "a) Using PSPG stabilization yields a conservative formulation (Hughes & Wells 2005)\n"
             "b) Instablities may occur for complex flow situations"
          << std::endl;
    }
  }

  // set flag for physical type of fluid flow
  physicaltype_ = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params, "Physical Type");
  if (((physicaltype_ == INPAR::FLUID::loma) or
          (physicaltype_ == INPAR::FLUID::varying_density)) and
      (fldparatimint_->IsStationary() == true))
    dserror("physical type is not supported in stationary FLUID implementation.");

  // set flag for type of linearization (fixed-point-like or Newton)
  //  fix-point like for Oseen or Stokes problems
  if (DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(params, "Linearisation") ==
      INPAR::FLUID::Newton)
  {
    if ((physicaltype_ == INPAR::FLUID::oseen) or (physicaltype_ == INPAR::FLUID::stokes))
      dserror(
          "Full Newton-linearization does not make sense for Oseen or Stokes problems.\nThey are "
          "already linear problems. Fix input file!");
    is_newton_ = true;
  }

  if (fldparatimint_->IsGenalphaNP() and physicaltype_ == INPAR::FLUID::loma)
    dserror("the combination Np_Gen_Alpha and loma is not supported");

  if (not fldparatimint_->IsGenalpha() and physicaltype_ == INPAR::FLUID::loma)
    dserror("the combination OST and loma is said to be supported but does not work!!");

  if (fldparatimint_->IsGenalphaNP() and is_conservative_)
    dserror("the combination Np_Gen_Alpha and conservative flow is not supported");

  if (not fldparatimint_->IsStationary() and is_conservative_ and
      physicaltype_ != INPAR::FLUID::incompressible)
  {
    if (myrank == 0)
      std::cout << std::endl
                << "Warning: missing time derivative terms in conservative formulation for "
                   "variable density flows!"
                << std::endl;
  }

  // get function number of given Oseen advective field if necessary
  if (physicaltype_ == INPAR::FLUID::oseen)
    oseenfieldfuncno_ = DRT::INPUT::get<int>(params, "OSEENFIELDFUNCNO");

  // ---------------------------------------------------------------------
  // get control parameters for stabilization and higher-order elements
  //----------------------------------------------------------------------
  Teuchos::ParameterList& stablist = params.sublist("RESIDUAL-BASED STABILIZATION");

  stabtype_ = DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(stablist, "STABTYPE");

  if (stabtype_ == INPAR::FLUID::stabtype_residualbased)
  {
    // no safety check necessary since all options are used
    tds_ = DRT::INPUT::IntegralValue<INPAR::FLUID::SubscalesTD>(stablist, "TDS");
    transient_ = DRT::INPUT::IntegralValue<INPAR::FLUID::Transient>(stablist, "TRANSIENT");
    pspg_ = DRT::INPUT::IntegralValue<int>(stablist, "PSPG");
    supg_ = DRT::INPUT::IntegralValue<int>(stablist, "SUPG");
    vstab_ = DRT::INPUT::IntegralValue<INPAR::FLUID::VStab>(stablist, "VSTAB");
    rstab_ = DRT::INPUT::IntegralValue<INPAR::FLUID::RStab>(stablist, "RSTAB");
    graddiv_ = DRT::INPUT::IntegralValue<int>(stablist, "GRAD_DIV");
    cross_ = DRT::INPUT::IntegralValue<INPAR::FLUID::CrossStress>(stablist, "CROSS-STRESS");
    reynolds_ =
        DRT::INPUT::IntegralValue<INPAR::FLUID::ReynoldsStress>(stablist, "REYNOLDS-STRESS");

    if (supg_ and (physicaltype_ == INPAR::FLUID::stokes))
      dserror(
          "Having SUPG-stabilization switched on (by default?) for Stokes problems, does not make "
          "sense! Please turn on brain before using BACI!");

    // overrule higher_order_ele if input-parameter is set
    // this might be interesting for fast (but slightly
    // less accurate) computations
    is_inconsistent_ = DRT::INPUT::IntegralValue<int>(stablist, "INCONSISTENT");

    is_reconstructder_ = DRT::INPUT::IntegralValue<int>(stablist, "Reconstruct_Sec_Der");
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
            whichtau_ == INPAR::FLUID::tau_codina or
            whichtau_ == INPAR::FLUID::tau_codina_convscaled or
            whichtau_ == INPAR::FLUID::tau_codina_wo_dt or
            whichtau_ == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina or
            whichtau_ == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt or
            whichtau_ == INPAR::FLUID::tau_hughes_franca_balestra_wo_dt))
      dserror("Definition of Tau cannot be handled by the element");

    // set correct stationary definition of stabilization parameter automatically
    if (fldparatimint_->IsStationary())
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
      else if (whichtau_ == INPAR::FLUID::tau_codina_convscaled)
        dserror("This stabilization Parameter is not implemented for stationary flows");
    }

    // tau_donea_huerta_wo_dt should only be used for stokes problems and vice versa
    if ((physicaltype_ == INPAR::FLUID::stokes) !=
        (whichtau_ == INPAR::FLUID::tau_hughes_franca_balestra_wo_dt))
      dserror(
          "Tau from Hughes, Franca & Balestra should only be used for Stokes problems and vice "
          "versa!");

    // get and check characteristic element length for stabilization parameter tau_Mu
    charelelengthu_ =
        DRT::INPUT::IntegralValue<INPAR::FLUID::CharEleLengthU>(stablist, "CHARELELENGTH_U");
    if (not(charelelengthu_ == INPAR::FLUID::streamlength_u or
            charelelengthu_ == INPAR::FLUID::volume_equivalent_diameter_u or
            charelelengthu_ == INPAR::FLUID::root_of_volume_u))
      dserror("Unknown characteristic element length for tau_Mu!");

    // get and check characteristic element length for stabilization parameter
    // tau_Mp and tau_C
    charelelengthpc_ =
        DRT::INPUT::IntegralValue<INPAR::FLUID::CharEleLengthPC>(stablist, "CHARELELENGTH_PC");
    if (not(charelelengthpc_ == INPAR::FLUID::streamlength_pc or
            charelelengthpc_ == INPAR::FLUID::volume_equivalent_diameter_pc or
            charelelengthpc_ == INPAR::FLUID::root_of_volume_pc))
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
      if (rstab_ == INPAR::FLUID::reactive_stab_usfem)
        viscreastabfac_ = -1.0;
      else if (rstab_ == INPAR::FLUID::reactive_stab_gls)
        viscreastabfac_ = 1.0;
    }

    // XFEM specific ghost penalty stabilization are set in the fluid_ele_parameter_intface-class
  }
  else if (stabtype_ == INPAR::FLUID::stabtype_edgebased)
  {
    if (myrank == 0)
    {
      IO::cout << "+-------------------------------------------------------------------------------"
                  "---+\n";
      IO::cout << " Edge-based stabilization: all residual-based stabilization terms are switched "
                  "off!\n";
      IO::cout << "+-------------------------------------------------------------------------------"
                  "---+\n"
               << IO::endl;
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

    // all edge-based flags are set in the fluid_ele_parameter_intface-class
  }
  else if (stabtype_ == INPAR::FLUID::stabtype_pressureprojection)
  {
    if (not(fldparatimint_->IsStationary() and
            ((physicaltype_ == INPAR::FLUID::stokes) or (physicaltype_ == INPAR::FLUID::oseen))))
      dserror(
          "Polynomial pressure projection has only been tested for stationary Stokes/Oseen "
          "problems. \n"
          "But it should work for other problems as well but only to circumvent "
          "inf-sup-instabilities! \n"
          "Convection instabilities have to be accounted for. \n"
          "Note that for now all residual-based stabilizations are switched off. \n"
          "Remove this dserror at own risk and have fun!");

    if (myrank == 0)
    {
      IO::cout << "+-------------------------------------------------------------------------------"
                  "---+\n";
      IO::cout << " Polynomial pressure projection: no residual-based, in particular no convective "
                  "stabilization! \n";
      IO::cout << "+-------------------------------------------------------------------------------"
                  "---+\n"
               << IO::endl;
    }
    //---------------------------------
    // if polynomial pressure projection stabilization is selected, all
    // residual-based stabilization terms except for SUGP are switched off


    pspg_ = false;
    supg_ = false;
    graddiv_ = false;
    vstab_ = INPAR::FLUID::viscous_stab_none;
    rstab_ = INPAR::FLUID::reactive_stab_none;
    cross_ = INPAR::FLUID::cross_stress_stab_none;
    reynolds_ = INPAR::FLUID::reynolds_stress_stab_none;
    tds_ = INPAR::FLUID::subscales_quasistatic;
    transient_ = INPAR::FLUID::inertia_stab_drop;
    is_inconsistent_ = false;

    //---------------------------------
    // polynomial pressure projection is switched on
    ppp_ = true;
  }
  else if (stabtype_ == INPAR::FLUID::stabtype_nostab)
  {
    if (myrank == 0)
    {
      IO::cout << "+-------------------------------------------------------------------------------"
                  "---+\n";
      IO::cout << "+                                   WARNING\n";
      IO::cout << " No stabilization selected: all stabilization terms are switched off!\n";
      IO::cout << "                            BACI says: Good luck!\n";
      IO::cout << "+-------------------------------------------------------------------------------"
                  "---+\n"
               << IO::endl;
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
  if ((tds_ == INPAR::FLUID::subscales_time_dependent) or
      (transient_ != INPAR::FLUID::inertia_stab_drop))
  {
    if (not fldparatimint_->IsGenalphaNP())
      dserror(
          "time dependent subscales does not work for OST/AfGenAlpha/BDF2/Stationary. \nOne need "
          "to look for bugs");
  }

  //---------------------------------
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

  return;
}

//----------------------------------------------------------------------*
//  set loma parameters                                  rasthofer 03/12|
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementLomaParameter(Teuchos::ParameterList& params)
{
  // get parameter lists
  Teuchos::ParameterList& lomaparams = params.sublist("LOMA");
  Teuchos::ParameterList& stabparams = params.sublist("RESIDUAL-BASED STABILIZATION");
  Teuchos::ParameterList& turbmodelparamsmfs = params.sublist("MULTIFRACTAL SUBGRID SCALES");

  //---------------------------------------------------------------------------------
  // material update with subgrid-scale temperature
  //---------------------------------------------------------------------------------

  update_mat_ = lomaparams.get<bool>("update material", false);

  //---------------------------------------------------------------------------------
  // parameter for additional rbvmm terms in continuity equation
  //---------------------------------------------------------------------------------

  conti_supg_ = DRT::INPUT::IntegralValue<int>(stabparams, "LOMA_CONTI_SUPG");
  conti_cross_ =
      DRT::INPUT::IntegralValue<INPAR::FLUID::CrossStress>(stabparams, "LOMA_CONTI_CROSS_STRESS");
  conti_reynolds_ = DRT::INPUT::IntegralValue<INPAR::FLUID::ReynoldsStress>(
      stabparams, "LOMA_CONTI_REYNOLDS_STRESS");

  //---------------------------------------------------------------------------------
  // parameter for additional multifractal subgrid-scale terms
  //---------------------------------------------------------------------------------

  if (turb_mod_action_ == INPAR::FLUID::multifractal_subgrid_scales)
    multifrac_loma_conti_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs, "LOMA_CONTI");

  return;
}

//----------------------------------------------------------------------*
//  set two phase parameters                                winter 05/14|
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementTwoPhaseParameter(Teuchos::ParameterList& params)
{
  // Two Phase Flow specific parameters,
  // Smeared specific parameters
  Teuchos::ParameterList& smearedlist = params.sublist("SMEARED");
  interface_thickness_ = smearedlist.get<double>("INTERFACE_THICKNESS");
  enhanced_gaussrule_ = DRT::INPUT::IntegralValue<int>(smearedlist, "ENHANCED_GAUSSRULE");

  // Surface tension specific parameters
  Teuchos::ParameterList& surftenslist = params.sublist("SURFACE TENSION");
  eval_surfacetension_ = DRT::INPUT::IntegralValue<INPAR::TWOPHASE::SurfaceTensionApprox>(
      surftenslist, "SURFTENSAPPROX");

  if (eval_surfacetension_ != INPAR::TWOPHASE::surface_tension_approx_none)
    include_surface_tension_ = true;

  return;
}


//----------------------------------------------------------------------*
//  set topopt parameters                               winklmaier 07/13|
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementTopoptParameter(Teuchos::ParameterList& params)
{
  topopt_params_[0] = params.get<double>("MIN_PORO");
  topopt_params_[1] = params.get<double>("MAX_PORO");
  topopt_params_[2] = params.get<double>("SMEAR_FAC");
  reaction_ = true;
  reaction_topopt_ = true;

  return;
}

//----------------------------------------------------------------------*
//  set turbulence parameters                            rasthofer 11/11|
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementTurbulenceParameters(
    Teuchos::ParameterList& params)
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
    const std::string fssgvdef = turbmodelparams.get<std::string>("FSSUGRVISC", "No");

    if (fssgvdef == "Smagorinsky_all")
      fssgv_ = INPAR::FLUID::smagorinsky_all;
    else if (fssgvdef == "Smagorinsky_small")
      fssgv_ = INPAR::FLUID::smagorinsky_small;
  }

  // get Smagorinsky model parameter for fine-scale subgrid viscosity
  // (Since either all-scale Smagorinsky model (i.e., classical LES model
  // as will be inititalized below) or fine-scale Smagorinsky model is
  // used (and never both), the same input parameter can be exploited.)
  if (fssgv_ != INPAR::FLUID::no_fssgv)
    Cs_ = turbmodelparamssgvisc.get<double>("C_SMAGORINSKY", 0.0);

  //---------------------------------------------------------------------------------
  // parameter for turbulence approach
  //---------------------------------------------------------------------------------

  // the default action is no model
  turb_mod_action_ = INPAR::FLUID::no_model;

  // No turbulent flow: TURBULENCE_APPROACH = DNS
  if (turbmodelparams.get<std::string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
  {
    if (fldparatimint_->IsStationary() == true)
      dserror("Stationary turbulent flow does not make any sense");

    std::string& physical_turbulence_model = turbmodelparams.get<std::string>("PHYSICAL_MODEL");

    // --------------------------------------------------
    // standard constant coefficient Smagorinsky model
    if (physical_turbulence_model == "Smagorinsky")
    {
      // the classic Smagorinsky model only requires one constant parameter
      turb_mod_action_ = INPAR::FLUID::smagorinsky;
      Cs_ = turbmodelparamssgvisc.get<double>("C_SMAGORINSKY");
      include_Ci_ = DRT::INPUT::IntegralValue<int>(turbmodelparamssgvisc, "C_INCLUDE_CI");
      Ci_ = turbmodelparamssgvisc.get<double>("C_YOSHIZAWA");
    }
    // --------------------------------------------------
    // Smagorinsky model with van Driest damping
    else if (physical_turbulence_model == "Smagorinsky_with_van_Driest_damping")
    {
      // that's only implemented for turbulent channel flow
      // wall function length is hard implemented
      if (turbmodelparamssgvisc.get<std::string>("CANONICAL_FLOW", "no") !=
          "channel_flow_of_height_2")
        dserror("van_Driest_damping only for channel_flow_of_height_2\n");

      // for the Smagorinsky model with van Driest damping, we need
      // a viscous length to determine the y+ (heigth in wall units)
      turb_mod_action_ = INPAR::FLUID::smagorinsky_with_van_Driest_damping;

      // get parameters of model
      Cs_ = turbmodelparamssgvisc.get<double>("C_SMAGORINSKY");
      l_tau_ = turbmodelparamssgvisc.get<double>("CHANNEL_L_TAU");
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
      Cs_averaged_ =
          DRT::INPUT::IntegralValue<int>(turbmodelparamssgvisc, "C_SMAGORINSKY_AVERAGED");
      Ci_ = turbmodelparamssgvisc.get<double>("C_YOSHIZAWA");
      include_Ci_ = DRT::INPUT::IntegralValue<int>(turbmodelparamssgvisc, "C_INCLUDE_CI");
    }
    else if (physical_turbulence_model == "Multifractal_Subgrid_Scales")
    {
      turb_mod_action_ = INPAR::FLUID::multifractal_subgrid_scales;

      // get parameters of model
      Csgs_ = turbmodelparamsmfs.get<double>("CSGS");
      Csgs_phi_ = turbmodelparamsmfs.get<double>("CSGS_PHI");
      adapt_Csgs_phi_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs, "ADAPT_CSGS_PHI");

      if (turbmodelparamsmfs.get<std::string>("SCALE_SEPARATION") == "algebraic_multigrid_operator")
        alpha_ = 3.0;
      else if (turbmodelparamsmfs.get<std::string>("SCALE_SEPARATION") == "box_filter")
        alpha_ = 2.0;
      else
        dserror("Unknown filter type!");

      CalcN_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs, "CALC_N");

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
      c_diff_ = turbmodelparamsmfs.get<double>("C_DIFF");  // loma only

      near_wall_limit_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs, "NEAR_WALL_LIMIT");
      near_wall_limit_scatra_ =
          DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs, "NEAR_WALL_LIMIT_CSGS_PHI");

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

      consistent_mfs_residual_ =
          DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs, "CONSISTENT_FLUID_RESIDUAL");
    }
    else if (physical_turbulence_model == "Vreman")
    {
      turb_mod_action_ = INPAR::FLUID::vreman;
      Cs_ = turbmodelparamssgvisc.get<double>("C_SMAGORINSKY");

      if (turbmodelparamssgvisc.get<std::string>("FILTER_WIDTH", "CubeRootVol") ==
          "Direction_dependent")
        vrfi_ = INPAR::FLUID::dir_dep;
      else if (turbmodelparamssgvisc.get<std::string>("FILTER_WIDTH", "CubeRootVol") ==
               "Minimum_length")
        vrfi_ = INPAR::FLUID::min_len;
      else
        vrfi_ = INPAR::FLUID::cuberootvol;
    }
    else if (physical_turbulence_model == "Dynamic_Vreman")
    {
      turb_mod_action_ = INPAR::FLUID::dynamic_vreman;
      if (turbmodelparamssgvisc.get<std::string>("FILTER_WIDTH", "CubeRootVol") ==
          "Direction_dependent")
        vrfi_ = INPAR::FLUID::dir_dep;
      else if (turbmodelparamssgvisc.get<std::string>("FILTER_WIDTH", "CubeRootVol") ==
               "Minimum_length")
        vrfi_ = INPAR::FLUID::min_len;
      else
        vrfi_ = INPAR::FLUID::cuberootvol;
    }
    else
    {
      dserror("Up to now, only Smagorinsky, Vreman and Multifractal Subgrid Scales are available");
    }
  }  // end if(Classical LES)
}
