/*----------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static parameters required for the evaluation of a standard
scalar transport element

This singleton class holds all static parameters required for the evaluation of a standard scalar
transport element, e.g., stabilization parameters and finite difference check parameters. All
parameters are usually set only once at the beginning of a simulation, namely during initialization
of the global time integrator, and then never touched again throughout the simulation. Enhanced
scalar transport problems, such as electrochemistry and levelset problems, instantiate additional,
problem specific singleton classes holding additional static parameters required for element
evaluation. These additional singleton classes are not meant to be derived from, but rather to
coexist with this general class.

\level 1

*/
/*----------------------------------------------------------------------*/
#include "baci_scatra_ele_parameter_std.H"

#include "baci_scatra_ele_parameter_timint.H"
#include "baci_utils_exceptions.H"
#include "baci_utils_singleton_owner.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterStd* DRT::ELEMENTS::ScaTraEleParameterStd::Instance(
    const std::string& disname  //!< name of discretization
)
{
  static auto singleton_map =
      CORE::UTILS::MakeSingletonMap<std::string>([](const std::string& disname)
          { return std::unique_ptr<ScaTraEleParameterStd>(new ScaTraEleParameterStd(disname)); });

  return singleton_map[disname].Instance(CORE::UTILS::SingletonAction::create, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterStd::ScaTraEleParameterStd(
    const std::string& disname  //!< name of discretization
    )
    : is_ale_(false),
      is_conservative_(false),
      sphericalcoords_(false),
      calcflux_domain_(INPAR::SCATRA::flux_none),
      writefluxids_(Teuchos::null),
      fdcheck_(INPAR::SCATRA::fdcheck_none),
      fdcheckeps_(0.),
      fdchecktol_(0.),
      nds_disp_(-1),
      nds_growth_(-1),
      nds_micro_(-1),
      nds_pres_(-1),
      nds_scatra_(-1),
      nds_thermo_(-1),
      nds_two_tensor_quantitiy_(-1),
      nds_vel_(-1),
      nds_wss_(-1),
      probnum_(0),
      semiimplicit_(false),
      intlayergrowth_convtol_(0.),
      intlayergrowth_itemax_(0),
      partitioned_multiscale_(false),
      is_emd_(false),
      emd_source_(-1),
      stabtype_(INPAR::SCATRA::stabtype_no_stabilization),
      whichtau_(INPAR::SCATRA::tau_zero),
      charelelength_(INPAR::SCATRA::streamlength),
      diffreastafac_(0.0),
      sgvel_(false),
      assgd_(false),
      whichassgd_(INPAR::SCATRA::assgd_artificial),
      tau_gp_(false),
      mat_gp_(false),
      tau_value_(0.),
      // we have to know the time parameters here to check for illegal combinations
      scatraparatimint_(DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(disname))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterStd::SetNodesetParameters(Teuchos::ParameterList& parameters)
{
  nds_disp_ = parameters.get<int>("ndsdisp", -1);
  nds_growth_ = parameters.get<int>("ndsgrowth", -1);
  nds_pres_ = parameters.get<int>("ndspres", -1);
  nds_scatra_ = parameters.get<int>("ndsscatra", -1);
  nds_thermo_ = parameters.get<int>("ndsthermo", -1);
  nds_two_tensor_quantitiy_ = parameters.get<int>("ndsTwoTensorQuantity", -1);
  nds_vel_ = parameters.get<int>("ndsvel", -1);
  nds_wss_ = parameters.get<int>("ndswss", -1);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterStd::SetParameters(Teuchos::ParameterList& parameters)
{
  // set ale case
  is_ale_ = parameters.get<bool>("isale", false);

  // set flag for conservative form
  const INPAR::SCATRA::ConvForm convform =
      DRT::INPUT::get<INPAR::SCATRA::ConvForm>(parameters, "convform");

  is_conservative_ = false;
  if (convform == INPAR::SCATRA::convform_conservative) is_conservative_ = true;

  // flag for writing the flux vector fields
  calcflux_domain_ = DRT::INPUT::get<INPAR::SCATRA::FluxType>(parameters, "calcflux_domain");

  //! vector containing ids of scalars for which flux vectors are calculated
  if (calcflux_domain_ != INPAR::SCATRA::flux_none)
    writefluxids_ = parameters.get<Teuchos::RCP<std::vector<int>>>("writeflux_ids");

  // set parameters for stabilization
  Teuchos::ParameterList& stablist = parameters.sublist("stabilization");

  // get definition for stabilization parameter tau
  whichtau_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(stablist, "DEFINITION_TAU");

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

  if (whichtau_ == INPAR::SCATRA::tau_numerical_value)
    tau_value_ = parameters.sublist("stabilization").get<double>("TAU_VALUE");

  // get characteristic element length for stabilization parameter definition
  charelelength_ =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::CharEleLength>(stablist, "CHARELELENGTH");

  // set (sign) factor for diffusive and reactive stabilization terms
  // (factor is zero for SUPG) and overwrite tau definition when there
  // is no stabilization
  stabtype_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(stablist, "STABTYPE");
  switch (stabtype_)
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
    case INPAR::SCATRA::stabtype_hdg_centered:
    case INPAR::SCATRA::stabtype_hdg_upwind:
      if (whichtau_ != INPAR::SCATRA::tau_numerical_value or tau_value_ <= 0.0)
        dserror(
            "Wrong definition for tau for hdg stabilization, only tau_numerical_value is allowed "
            "with tau>0");
      break;
    default:
      dserror("unknown definition for stabilization parameter");
  }

  // set flags for subgrid-scale velocity and all-scale subgrid-diffusivity term
  // (default: "false" for both flags)
  sgvel_ = DRT::INPUT::IntegralValue<int>(stablist, "SUGRVEL");
  assgd_ = DRT::INPUT::IntegralValue<int>(stablist, "ASSUGRDIFF");

  // select type of all-scale subgrid diffusivity if included
  whichassgd_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::AssgdType>(stablist, "DEFINITION_ASSGD");

  // set flags for potential evaluation of tau and material law at int. point
  const INPAR::SCATRA::EvalTau tauloc =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalTau>(stablist, "EVALUATION_TAU");
  tau_gp_ = (tauloc == INPAR::SCATRA::evaltau_integration_point);  // set true/false

  const INPAR::SCATRA::EvalMat matloc =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalMat>(stablist, "EVALUATION_MAT");
  mat_gp_ = (matloc == INPAR::SCATRA::evalmat_integration_point);  // set true/false

  // check for illegal combinations
  if (sgvel_ or assgd_)
  {
    // check for matching flags
    if (not mat_gp_ or not tau_gp_)
      dserror(
          "Evaluation of material and stabilization parameters need to be done at the integration "
          "points if subgrid-scale velocity is included!");
  }

  // get quantities for finite difference check
  fdcheck_ = DRT::INPUT::get<INPAR::SCATRA::FDCheck>(parameters, "fdcheck");
  fdcheckeps_ = parameters.get<double>("fdcheckeps");
  fdchecktol_ = parameters.get<double>("fdchecktol");

  // set flag for use of spherical coordinates
  sphericalcoords_ = parameters.get<bool>("sphericalcoords", false);

  // set problem number
  probnum_ = parameters.get<int>("probnum", 0);

  // set evaluation type
  semiimplicit_ = parameters.get<bool>("semiimplicit", false);

  // set local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving
  // interface layer growth
  intlayergrowth_convtol_ = parameters.get<double>("intlayergrowth_convtol", 0.);

  // set maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling
  // involving interface layer growth
  intlayergrowth_itemax_ = parameters.get<unsigned>("intlayergrowth_itemax", 0);

  // set flag for truly partitioned multi-scale simulation
  partitioned_multiscale_ = parameters.get<bool>("partitioned_multiscale", false);

  // set flag for electromagnetic diffusion simulation
  is_emd_ = parameters.get<bool>("electromagnetic_diffusion", false);

  // electromagnetic diffusion source function
  emd_source_ = parameters.get<int>("electromagnetic_diffusion_source", -1);
}

int DRT::ELEMENTS::ScaTraEleParameterStd::NdsDisp() const
{
  dsassert(nds_disp_ != -1,
      "You try to access the number of dofset associated with displacement dofs without "
      "having set it!");
  return nds_disp_;
}

int DRT::ELEMENTS::ScaTraEleParameterStd::NdsGrowth() const
{
  dsassert(nds_growth_ != -1,
      "You try to access the number of dofset associated with interface growth dofs without "
      "having set it!");
  return nds_growth_;
}

int DRT::ELEMENTS::ScaTraEleParameterStd::NdsMicro() const
{
  dsassert(nds_micro_ != -1,
      "You try to access the number of dofset to write micro scale values on without having "
      "set it!");
  return nds_micro_;
}

int DRT::ELEMENTS::ScaTraEleParameterStd::NdsPres() const
{
  dsassert(nds_pres_ != -1,
      "You try to access the number of dofset associated with pressure dofs without having "
      "set it!");
  return nds_pres_;
}

int DRT::ELEMENTS::ScaTraEleParameterStd::NdsScaTra() const
{
  dsassert(nds_scatra_ != -1,
      "You try to access the number of dofset associated with scalar transport dofs without having "
      "set it!");
  return nds_scatra_;
}

int DRT::ELEMENTS::ScaTraEleParameterStd::NdsThermo() const
{
  dsassert(nds_thermo_ != -1,
      "You try to access the number of dofset associated with temperature dofs without having "
      "set it!");
  return nds_thermo_;
}

int DRT::ELEMENTS::ScaTraEleParameterStd::NdsTwoTensorQuantity() const
{
  dsassert(nds_two_tensor_quantitiy_ != -1,
      "You try to access the number of dofset associated with two-tensor quantity dofs without "
      "having set it!");
  return nds_two_tensor_quantitiy_;
}

int DRT::ELEMENTS::ScaTraEleParameterStd::NdsVel() const
{
  dsassert(nds_vel_ != -1,
      "You try to access the number of dofset associated with velocity related dofs without "
      "having set it!");
  return nds_vel_;
}

int DRT::ELEMENTS::ScaTraEleParameterStd::NdsWss() const
{
  dsassert(nds_wss_ != -1,
      "You try to access the number of dofset associated with wall shear stress dofs without "
      "having set it!");
  return nds_wss_;
}
BACI_NAMESPACE_CLOSE
