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
  is_ale_(false),
  is_conservative_(false),
  writeflux_(INPAR::SCATRA::flux_no),
  writefluxids_(Teuchos::null),
  fdcheck_(INPAR::SCATRA::fdcheck_none),
  fdcheckeps_(0.),
  fdchecktol_(0.),
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
  scalarforcing_(INPAR::FLUID::scalarforcing_no),
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
{
  // we have to know the time parameters here to check for illegal combinations
  scatraparatimint_ = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance();

  return;
}


//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameter::SetElementGeneralParameters(Teuchos::ParameterList& params)
{
  // set ale case
  is_ale_ = params.get<bool>("isale",false);

  // set flag for conservative form
  const INPAR::SCATRA::ConvForm convform =
    DRT::INPUT::get<INPAR::SCATRA::ConvForm>(params, "form of convective term");
  is_conservative_ = false;
  if (convform ==INPAR::SCATRA::convform_conservative) is_conservative_ = true;

  // flag for writing the flux vector fields
  writeflux_ =  DRT::INPUT::get<INPAR::SCATRA::FluxType>(params, "writeflux");

  //! vector containing ids of scalars for which flux vectors are calculated
  writefluxids_ =  params.get<Teuchos::RCP<std::vector<int> > >("writefluxids");

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

  // get quantities for finite difference check
  fdcheck_ = DRT::INPUT::get<INPAR::SCATRA::FDCheck>(params,"fdcheck");
  fdcheckeps_ = params.get<double>("fdcheckeps");
  fdchecktol_ = params.get<double>("fdchecktol");

  return;
}


//----------------------------------------------------------------------*
//  set turbulence parameters                            rasthofer 11/11|
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameter::SetElementTurbulenceParameters( Teuchos::ParameterList& params )
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

  // define forcing for scalar field
  if (turbulencelist.get<std::string>("SCALAR_FORCING","no")=="isotropic")
    scalarforcing_ = INPAR::FLUID::scalarforcing_isotropic;
  else if(turbulencelist.get<std::string>("SCALAR_FORCING","no")=="mean_scalar_gradient")
    scalarforcing_ = INPAR::FLUID::scalarforcing_mean_scalar_gradient;
  else
    scalarforcing_ = INPAR::FLUID::scalarforcing_no;

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

  // in some cases we may want to switch off the turbulence model in the scalar field
  if (not DRT::INPUT::IntegralValue<int>(turbulencelist,"TURBMODEL_LS"))
  {
    fssgd_ = false;
    whichfssgd_ = INPAR::SCATRA::fssugrdiff_no;
    turbmodel_ = INPAR::FLUID::no_model;
  }

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

      adapt_Csgs_phi_ = DRT::INPUT::IntegralValue<bool>(mfslist,"ADAPT_CSGS_PHI");

      // safety  check
      if (BD_gp_ and (not mat_gp_)) dserror("evaluation of B and D at gauss-point should always be combined with evaluation material at gauss-point!");
    }
  }

  turbinflow_ = params.get<bool>("turbulent inflow");

  return;
}

