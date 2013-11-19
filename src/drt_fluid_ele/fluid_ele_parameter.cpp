/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_std.cpp

\brief Evaluation of general fluid parameter for standard fluid

       Poro specific parameters are defined in a derived class.
       Since there are a couple of terms from the stdfluid called
       by the porofluid this is  inevitable.

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/
#include "fluid_ele_parameter.H"
#include "fluid_ele_parameter_base.H"
#include <string>
#include <iostream>
#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io_pstream.H"

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameter::FluidEleParameter()
  : DRT::ELEMENTS::FluidEleParameterBase::FluidEleParameterBase(),
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
    tau_gp_(false),      // standard evaluation of tau at the element center
    turb_mod_action_(INPAR::FLUID::no_model), // turbulence parameters
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
    meanCai_(0.0)
{

}

//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameter::SetElementStdFluidParameter( Teuchos::ParameterList& params,
                                                                        int myrank )
{
  // call method from base class
  SetElementGeneralFluidParameter(params,myrank);

 if (fldparatimint_->IsGenalphaNP() and physicaltype_ == INPAR::FLUID::loma)
   dserror("the combination Np_Gen_Alpha and loma is not supported");

 if (not fldparatimint_->IsGenalpha() and physicaltype_ == INPAR::FLUID::loma)
   dserror("the combination OST and loma is said to be supported but does not work!!");

 if (fldparatimint_->IsGenalphaNP() and is_conservative_)
   dserror("the combination Np_Gen_Alpha and conservative flow is not supported");

 if (not fldparatimint_->IsStationary() and is_conservative_ and physicaltype_ != INPAR::FLUID::incompressible)
 {
   if (myrank == 0)
     std::cout << std::endl << "Warning: missing time derivative terms in conservative formulation for variable density flows!" << std::endl;
 }

 // ---------------------------------------------------------------------
 // get control parameters for stabilization and higher-order elements
 //----------------------------------------------------------------------
 Teuchos::ParameterList& stablist = params.sublist("RESIDUAL-BASED STABILIZATION");
 Teuchos::ParameterList& stablist_edgebased = params.sublist("EDGE-BASED STABILIZATION");

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
   graddiv_  = DRT::INPUT::IntegralValue<int>(stablist,"GRAD_DIV");
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
   if (not fldparatimint_->IsGenalphaNP())
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
    if (fldparatimint_->IsStationary() == true)
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
    else if (physical_turbulence_model == "Vreman")
    {
      turb_mod_action_ = INPAR::FLUID::vreman;
      Cs_ = turbmodelparamssgvisc.get<double>("C_SMAGORINSKY");

//      const std::string fssgvdef = turbmodelparams.get<std::string>("FSSUGRVISC","No");
//
//      if (fssgvdef == "Smagorinsky_all")        fssgv_ = INPAR::FLUID::smagorinsky_all;
//      else if (fssgvdef == "Smagorinsky_small") fssgv_ = INPAR::FLUID::smagorinsky_small;
      if(turbmodelparamssgvisc.get<std::string>("FILTER_WIDTH","CubeRootVol")=="Direction_dependent")
        vrfi_ = INPAR::FLUID::dir_dep;
      else if(turbmodelparamssgvisc.get<std::string>("FILTER_WIDTH","CubeRootVol")=="Minimum_length")
        vrfi_ = INPAR::FLUID::min_len;
      else
        vrfi_ = INPAR::FLUID::cuberootvol;

    }
        else if (physical_turbulence_model == "Dynamic_Vreman")
    {
      turb_mod_action_ = INPAR::FLUID::dynamic_vreman;
      if(turbmodelparamssgvisc.get<std::string>("FILTER_WIDTH","CubeRootVol")=="Direction_dependent")
        vrfi_ = INPAR::FLUID::dir_dep;
      else if(turbmodelparamssgvisc.get<std::string>("FILTER_WIDTH","CubeRootVol")=="Minimum_length")
        vrfi_ = INPAR::FLUID::min_len;
      else
        vrfi_ = INPAR::FLUID::cuberootvol;
    }
    else
    {
      dserror("Up to now, only Smagorinsky, Scale Similarity and Multifractal Subgrid Scales are available");
    }
  } // end if(Classical LES)
}

