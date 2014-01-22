/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter.cpp

\brief Setting of general scatra parameter for element evaluation

This file has to contain all parameters called in scatra_ele_calc.cpp.
Additional parameters required in derived classes of ScaTraEleCalc have to
be set in problem specific parameter lists derived from this class.

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_parameter.H"
//#include <string>
//#include <iostream>
#include "../drt_lib/drt_dserror.H"
//#include "../drt_io/io_pstream.H"


//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameter::ScaTraEleParameter()
 : //set_general_fluid_parameter_(false),
  set_general_scatra_parameter_(false),
  scatratype_(INPAR::SCATRA::scatratype_undefined),
  is_ale_(false),
  is_conservative_(false),
  stabtype_(INPAR::SCATRA::stabtype_no_stabilization),
  whichtau_(INPAR::SCATRA::tau_zero),
  charelelength_(INPAR::SCATRA::streamlength),
  diffreastafac_(0.0),
  sgvel_(false),
  assgd_(false),
  whichassgd_(INPAR::SCATRA::assgd_artificial),
  tau_gp_(false),
  mat_gp_(false),
  turbmodel_(INPAR::FLUID::no_model),
  fssgd_(false),
  whichfssgd_(INPAR::SCATRA::fssugrdiff_no),
  Cs_(0.0),
  tpn_(0.0),
  Cs_av_(false),
  Csgs_sgvel_(0.0),
  alpha_(0.0),
  calc_N_(false),
  N_vel_(0.0),
  refvel_(INPAR::FLUID::strainrate),
  reflength_(INPAR::FLUID::cube_edge),
  c_nu_(0.0),
  nwl_(false),
  nwl_scatra_(false),
  beta_(false),
  BD_gp_(false),
  Csgs_sgphi_(0.0),
  c_diff_(0.0),
  mfs_conservative_(false),
  meanCai_(0.0),
  adapt_Csgs_phi_(false),
  turbinflow_(false)
//    is_elch_((numdofpernode_ - numscal_) >= 1),  // bool set implicitely
//    is_reactive_(false),      // bool set
//    is_coupled_(false),            //bool set
//    diffreastafac_(0.0),       // set double (SUPG)
//    is_anisotropic_(false),   // bool set

//    sgvel_(false),            // bool set
//    betterconsistency_(false), // bool set
//    migrationintau_(true),     // bool set
//    migrationstab_(true),      // bool set
//    migrationinresidual_(true),// bool set
//    update_mat_(false),        // bool set
//    // whichtau_ not initialized
//    turbmodel_(INPAR::FLUID::no_model), // enum initialized
//    mfs_conservative_(false), // set false
//    HSTCConds_(0,NULL),   //vector containing homogeneous coupling conditions
//    epotnp_(true),      // initialized to zero
//    emagnetnp_(true),   // initialized to zero
//    gradpot_(true),     // initialized to zero
//    evelnp_(true),      // initialized to zero
//    econvelnp_(true),   // initialized to zero
//    efsvel_(true),      // initialized to zero
//    eaccnp_(true),      // initialized to zero
//    velint_(true),      // initialized to zero
//    convelint_(true),   // initialized to zero
//    sgvelint_(true),    // initialized to zero
//    fsvelint_(true),    // initialized to zero
//    mfsgvelint_(true),  // initialized to zero
//    migvelint_(true),   // initialized to zero
//    conv_(true),        // initialized to zero
//    sgconv_(true),      // initialized to zero
//    vdiv_(0.0),         // set double
//    mfsvdiv_(0.0),      // set double
//    eprenp_(true),      // initialized to zero
//    shc_(0.0),      // set double
//    visc_(0.0),     // set double
//    diffcond_(false),
//    cursolvar_(false),
//    chemdiffcoupltransp_(true),
//    chemdiffcouplcurr_(true),
//    constparams_(true),
//    newman_(false),
//    diffbased_(true),
//    gradphicoupling_(numscal_),
//    curdiv_(0.0),
//    trans_(numscal_,0.0),
//    transelim_(0.0),
//    transderiv_(numscal_,std::vector<double>(numscal_,0.0 )),
//    cond_(1,0.0),
//    condderiv_(numscal_,0.0),
//    diffusderiv_(numscal_,0.0),
//    diffuselimderiv_(numscal_,0.0),
//    diffuselim_(0.0),
//    eps_(1,1.0),
//    tort_(1,1.0),
//    epstort_(1,1.0),
//    a_(0.0),
//    b_(0.0),
//    c_(0.0),
//    ecurnp_(true),
//    curint_(true),
//    equpot_(INPAR::ELCH::equpot_enc),
//    frt_(0.0)
{
  // we have to know the time parameters here to check for illegal combinations
  scatraparatimint_ = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance();

  return;
}


//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameter::SetElementGeneralScaTraParameter(
  Teuchos::ParameterList& params,
  int myrank )
{
  //TODO: SCATRA_ELE_CLEANING: CalcInitTimeDeriv -> zweimal gerufen
  if(set_general_scatra_parameter_ == false)
    set_general_scatra_parameter_ = true;
  else
  {
    if (myrank == 0)
      std::cout << "General scatra parameters should be set only once!!\n " << std::endl;
  }


  // the type of scalar transport problem has to be provided for all actions!
 scatratype_ = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");
  if (scatratype_ == INPAR::SCATRA::scatratype_undefined)
    dserror("Set parameter SCATRATYPE in your input file!");

  // set ale case
  is_ale_ = params.get<bool>("isale",false);

  // set flag for conservative form
  const INPAR::SCATRA::ConvForm convform =
    DRT::INPUT::get<INPAR::SCATRA::ConvForm>(params, "form of convective term");
  is_conservative_ = false;
  if (convform ==INPAR::SCATRA::convform_conservative) is_conservative_ = true;

  // set parameters for stabilization
  Teuchos::ParameterList& stablist = params.sublist("STABILIZATION");

  // get definition for stabilization parameter tau
  whichtau_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(stablist,"DEFINITION_TAU");

  // set correct stationary definition for stabilization parameter automatically
  // and ensure that exact stabilization parameter is only used in stationary case
  if (scatraparatimint_->IsStationary())
  {
    if (whichtau_ == INPAR::SCATRA::tau_taylor_hughes_zarins)
      whichtau_ = INPAR::SCATRA::tau_taylor_hughes_zarins_wo_dt;
    else if (whichtau_ == INPAR::SCATRA::tau_franca_valentin)
      whichtau_ = INPAR::SCATRA::tau_franca_valentin_wo_dt;
    else if (whichtau_ == INPAR::SCATRA::tau_shakib_hughes_codina)
      whichtau_ = INPAR::SCATRA::tau_shakib_hughes_codina_wo_dt;
    else if (whichtau_ == INPAR::SCATRA::tau_codina)
      whichtau_ = INPAR::SCATRA::tau_codina_wo_dt;
    else if (whichtau_ == INPAR::SCATRA::tau_franca_madureira_valentin)
      whichtau_ = INPAR::SCATRA::tau_franca_madureira_valentin_wo_dt;
  }
  else
  {
    if (whichtau_ == INPAR::SCATRA::tau_exact_1d)
      dserror("exact stabilization parameter only available for stationary case");
  }

  // get characteristic element length for stabilization parameter definition
  charelelength_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::CharEleLength>(stablist,"CHARELELENGTH");

  // set (sign) factor for diffusive and reactive stabilization terms
  // (factor is zero for SUPG) and overwrite tau definition when there
  // is no stabilization
  stabtype_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(stablist,"STABTYPE");
  switch(stabtype_)
  {
  case INPAR::SCATRA::stabtype_no_stabilization:
    whichtau_ = INPAR::SCATRA::tau_zero;
    break;
  case INPAR::SCATRA::stabtype_SUPG:
    diffreastafac_ = 0.0;
    break;
  case INPAR::SCATRA::stabtype_GLS:
    diffreastafac_ = 1.0;
    break;
  case INPAR::SCATRA::stabtype_USFEM:
    diffreastafac_ = -1.0;
    break;
  default:
    dserror("unknown definition for stabilization parameter");
    break;
  }

  // set flags for subgrid-scale velocity and all-scale subgrid-diffusivity term
  // (default: "false" for both flags)
  sgvel_ = DRT::INPUT::IntegralValue<int>(stablist,"SUGRVEL");
  assgd_ = DRT::INPUT::IntegralValue<int>(stablist,"ASSUGRDIFF");

  // select type of all-scale subgrid diffusivity if included
  whichassgd_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::AssgdType>(stablist,"DEFINITION_ASSGD");

  // set flags for potential evaluation of tau and material law at int. point
  const INPAR::SCATRA::EvalTau tauloc = DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalTau>(stablist,"EVALUATION_TAU");
  tau_gp_ = (tauloc == INPAR::SCATRA::evaltau_integration_point); // set true/false
  const INPAR::SCATRA::EvalMat matloc = DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalMat>(stablist,"EVALUATION_MAT");
  mat_gp_ = (matloc == INPAR::SCATRA::evalmat_integration_point); // set true/false

  // check for illegal combinations
  if (sgvel_ or assgd_)
  {
    // check for matching flags
    if (not mat_gp_ or not tau_gp_)
     dserror("Evaluation of material and stabilization parameters need to be done at the integration points if subgrid-scale velocity is included!");
  }

  if (sgvel_ and scatratype_ == INPAR::SCATRA::scatratype_levelset)
    dserror("CalcSubgrVelocityLevelSet not available anymore");

  return;
}


//----------------------------------------------------------------------*
//  set turbulence parameters                            rasthofer 11/11|
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameter::SetElementTurbulenceParameter( Teuchos::ParameterList& params )
{
  // get list with model-specific parameters
  Teuchos::ParameterList& turbulencelist = params.sublist("TURBULENCE MODEL");
  Teuchos::ParameterList& sgvisclist = params.sublist("SUBGRID VISCOSITY");
  Teuchos::ParameterList& mfslist = params.sublist("MULTIFRACTAL SUBGRID SCALES");

  // set flag for turbulence model
  if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "no_model")
      turbmodel_ = INPAR::FLUID::no_model;
  else if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Smagorinsky")
    turbmodel_ = INPAR::FLUID::smagorinsky;
  else if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Dynamic_Smagorinsky")
    turbmodel_ = INPAR::FLUID::dynamic_smagorinsky;
  else if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Multifractal_Subgrid_Scales")
    turbmodel_ = INPAR::FLUID::multifractal_subgrid_scales;
  else if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Dynamic_Vreman")
    turbmodel_ = INPAR::FLUID::dynamic_vreman;
  else dserror("Unknown turbulence model for scatra!");

  // set flag for fine-scale subgrid diffusivity and perform some checks
  whichfssgd_ = DRT::INPUT::get<INPAR::SCATRA::FSSUGRDIFF>(params, "fs subgrid diffusivity");
  if (whichfssgd_ == INPAR::SCATRA::fssugrdiff_artificial)
  {
    fssgd_ = true;

    // check for solver type
    if (scatraparatimint_->IsIncremental()) dserror("Artificial fine-scale subgrid-diffusivity approach only in combination with non-incremental solver so far!");
  }
  else if (whichfssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_all or whichfssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
  {
    fssgd_ = true;

    // check for solver type
    if (not scatraparatimint_->IsIncremental()) dserror("Fine-scale subgrid-diffusivity approach using all/small-scale Smagorinsky model only in combination with incremental solver so far!");
  }

  // check for combination of all-scale and fine-scale subgrid diffusivity
  if (assgd_ and fssgd_) dserror("No combination of all-scale and fine-scale subgrid-diffusivity approach currently possible!");

  if (turbmodel_!=INPAR::FLUID::no_model or (scatraparatimint_->IsIncremental() and fssgd_))
  {
    // get Smagorinsky constant and turbulent Prandtl number
    Cs_  = sgvisclist.get<double>("C_SMAGORINSKY");
    tpn_ = sgvisclist.get<double>("C_TURBPRANDTL");
    if (tpn_ <= 1.0E-16) dserror("Turbulent Prandtl number should be larger than zero!");

    Cs_av_ = DRT::INPUT::IntegralValue<int>(sgvisclist,"C_SMAGORINSKY_AVERAGED");

    if (turbmodel_ == INPAR::FLUID::dynamic_vreman)
      Cs_ = turbulencelist.get<double>("Dt_vreman",1.0);

    // get model parameters
    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    {
      // necessary parameters for subgrid-scale velocity estimation
      Csgs_sgvel_ = mfslist.get<double>("CSGS");
      if (mfslist.get<std::string>("SCALE_SEPARATION") == "algebraic_multigrid_operator")
       alpha_ = 3.0;
      else dserror("Scale-Separtion method not supported!");
      calc_N_= DRT::INPUT::IntegralValue<int>(mfslist,"CALC_N");
      N_vel_ = mfslist.get<double>("N");
      if (mfslist.get<std::string>("REF_VELOCITY") == "strainrate")
       refvel_ = INPAR::FLUID::strainrate;
      else if (mfslist.get<std::string>("REF_VELOCITY") == "resolved")
       refvel_ = INPAR::FLUID::resolved;
      else if (mfslist.get<std::string>("REF_VELOCITY") == "fine_scale")
       refvel_ = INPAR::FLUID::fine_scale;
      else
       dserror("Unknown velocity!");
      if (mfslist.get<std::string>("REF_LENGTH") == "cube_edge")
       reflength_ = INPAR::FLUID::cube_edge;
      else if (mfslist.get<std::string>("REF_LENGTH") == "sphere_diameter")
       reflength_ = INPAR::FLUID::sphere_diameter;
      else if (mfslist.get<std::string>("REF_LENGTH") == "streamlength")
       reflength_ = INPAR::FLUID::streamlength;
      else if (mfslist.get<std::string>("REF_LENGTH") == "gradient_based")
       reflength_ = INPAR::FLUID::gradient_based;
      else if (mfslist.get<std::string>("REF_LENGTH") == "metric_tensor")
       reflength_ = INPAR::FLUID::metric_tensor;
      else
       dserror("Unknown length!");
      c_nu_ = mfslist.get<double>("C_NU");
      nwl_ = DRT::INPUT::IntegralValue<int>(mfslist,"NEAR_WALL_LIMIT");
      // necessary parameters for subgrid-scale scalar estimation
      Csgs_sgphi_ = mfslist.get<double>("CSGS_PHI");
      c_diff_ = mfslist.get<double>("C_DIFF");
      nwl_scatra_ = DRT::INPUT::IntegralValue<int>(mfslist,"NEAR_WALL_LIMIT_CSGS_PHI");
      // general parameters
      beta_ = mfslist.get<double>("BETA");
      if (beta_!=0.0) dserror("Lhs terms for mfs not included! Fixed-point iteration only!");
      if (mfslist.get<std::string>("EVALUATION_B") == "element_center")
      BD_gp_ = false;
      else if (mfslist.get<std::string>("EVALUATION_B") == "integration_point")
      BD_gp_ = true;
      else
        dserror("Unknown evaluation point!");
      if (mfslist.get<std::string>("CONVFORM") == "convective")
      mfs_conservative_ = false;
      else if (mfslist.get<std::string>("CONVFORM") == "conservative")
      mfs_conservative_ = true;
      else
        dserror("Unknown form of convective term!");
      if (scatratype_ == INPAR::SCATRA::scatratype_loma and mfs_conservative_)
        dserror("Conservative formulation not supported for loma!");

      adapt_Csgs_phi_ = DRT::INPUT::IntegralValue<bool>(mfslist,"ADAPT_CSGS_PHI");

      // safety  check
      if (BD_gp_ and (not mat_gp_)) dserror("evaluation of B and D at gauss-point should always be combined with evaluation material at gauss-point!");
    }
  }

  turbinflow_ = params.get<bool>("turbulent inflow");

  return;
}

