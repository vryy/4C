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
#include "4C_scatra_ele_parameter_std.hpp"

#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraEleParameterStd* Discret::ELEMENTS::ScaTraEleParameterStd::instance(
    const std::string& disname  //!< name of discretization
)
{
  static auto singleton_map =
      Core::UTILS::MakeSingletonMap<std::string>([](const std::string& disname)
          { return std::unique_ptr<ScaTraEleParameterStd>(new ScaTraEleParameterStd(disname)); });

  return singleton_map[disname].instance(Core::UTILS::SingletonAction::create, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraEleParameterStd::ScaTraEleParameterStd(
    const std::string& disname  //!< name of discretization
    )
    : is_ale_(false),
      is_conservative_(false),
      sphericalcoords_(false),
      calcflux_domain_(Inpar::ScaTra::flux_none),
      writefluxids_(Teuchos::null),
      fdcheck_(Inpar::ScaTra::fdcheck_none),
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
      has_external_force_(false),
      stabtype_(Inpar::ScaTra::stabtype_no_stabilization),
      whichtau_(Inpar::ScaTra::tau_zero),
      charelelength_(Inpar::ScaTra::streamlength),
      diffreastafac_(0.0),
      sgvel_(false),
      assgd_(false),
      whichassgd_(Inpar::ScaTra::assgd_artificial),
      tau_gp_(false),
      mat_gp_(false),
      tau_value_(0.),
      // we have to know the time parameters here to check for illegal combinations
      scatraparatimint_(Discret::ELEMENTS::ScaTraEleParameterTimInt::instance(disname))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterStd::set_nodeset_parameters(
    Teuchos::ParameterList& parameters)
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
void Discret::ELEMENTS::ScaTraEleParameterStd::set_parameters(Teuchos::ParameterList& parameters)
{
  // set ale case
  is_ale_ = parameters.get<bool>("isale", false);

  // set flag for conservative form
  const Inpar::ScaTra::ConvForm convform =
      Core::UTILS::GetAsEnum<Inpar::ScaTra::ConvForm>(parameters, "convform");

  is_conservative_ = false;
  if (convform == Inpar::ScaTra::convform_conservative) is_conservative_ = true;

  // flag for writing the flux vector fields
  calcflux_domain_ = Core::UTILS::GetAsEnum<Inpar::ScaTra::FluxType>(parameters, "calcflux_domain");

  //! vector containing ids of scalars for which flux vectors are calculated
  if (calcflux_domain_ != Inpar::ScaTra::flux_none)
    writefluxids_ = parameters.get<Teuchos::RCP<std::vector<int>>>("writeflux_ids");

  // set parameters for stabilization
  Teuchos::ParameterList& stablist = parameters.sublist("stabilization");

  // get definition for stabilization parameter tau
  whichtau_ = Core::UTILS::IntegralValue<Inpar::ScaTra::TauType>(stablist, "DEFINITION_TAU");

  // set correct stationary definition for stabilization parameter automatically
  // and ensure that exact stabilization parameter is only used in stationary case
  if (scatraparatimint_->is_stationary())
  {
    if (whichtau_ == Inpar::ScaTra::tau_taylor_hughes_zarins)
      whichtau_ = Inpar::ScaTra::tau_taylor_hughes_zarins_wo_dt;
    else if (whichtau_ == Inpar::ScaTra::tau_franca_valentin)
      whichtau_ = Inpar::ScaTra::tau_franca_valentin_wo_dt;
    else if (whichtau_ == Inpar::ScaTra::tau_shakib_hughes_codina)
      whichtau_ = Inpar::ScaTra::tau_shakib_hughes_codina_wo_dt;
    else if (whichtau_ == Inpar::ScaTra::tau_codina)
      whichtau_ = Inpar::ScaTra::tau_codina_wo_dt;
    else if (whichtau_ == Inpar::ScaTra::tau_franca_madureira_valentin)
      whichtau_ = Inpar::ScaTra::tau_franca_madureira_valentin_wo_dt;
  }
  else
  {
    if (whichtau_ == Inpar::ScaTra::tau_exact_1d)
      FOUR_C_THROW("exact stabilization parameter only available for stationary case");
  }

  if (whichtau_ == Inpar::ScaTra::tau_numerical_value)
    tau_value_ = parameters.sublist("stabilization").get<double>("TAU_VALUE");

  // get characteristic element length for stabilization parameter definition
  charelelength_ =
      Core::UTILS::IntegralValue<Inpar::ScaTra::CharEleLength>(stablist, "CHARELELENGTH");

  // set (sign) factor for diffusive and reactive stabilization terms
  // (factor is zero for SUPG) and overwrite tau definition when there
  // is no stabilization
  stabtype_ = Core::UTILS::IntegralValue<Inpar::ScaTra::StabType>(stablist, "STABTYPE");
  switch (stabtype_)
  {
    case Inpar::ScaTra::stabtype_no_stabilization:
      whichtau_ = Inpar::ScaTra::tau_zero;
      break;
    case Inpar::ScaTra::stabtype_SUPG:
      diffreastafac_ = 0.0;
      break;
    case Inpar::ScaTra::stabtype_GLS:
      diffreastafac_ = 1.0;
      break;
    case Inpar::ScaTra::stabtype_USFEM:
      diffreastafac_ = -1.0;
      break;
    case Inpar::ScaTra::stabtype_hdg_centered:
    case Inpar::ScaTra::stabtype_hdg_upwind:
      if (whichtau_ != Inpar::ScaTra::tau_numerical_value or tau_value_ <= 0.0)
        FOUR_C_THROW(
            "Wrong definition for tau for hdg stabilization, only tau_numerical_value is allowed "
            "with tau>0");
      break;
    default:
      FOUR_C_THROW("unknown definition for stabilization parameter");
  }

  // set flags for subgrid-scale velocity and all-scale subgrid-diffusivity term
  // (default: "false" for both flags)
  sgvel_ = Core::UTILS::IntegralValue<int>(stablist, "SUGRVEL");
  assgd_ = Core::UTILS::IntegralValue<int>(stablist, "ASSUGRDIFF");

  // select type of all-scale subgrid diffusivity if included
  whichassgd_ = Core::UTILS::IntegralValue<Inpar::ScaTra::AssgdType>(stablist, "DEFINITION_ASSGD");

  // set flags for potential evaluation of tau and material law at int. point
  const Inpar::ScaTra::EvalTau tauloc =
      Core::UTILS::IntegralValue<Inpar::ScaTra::EvalTau>(stablist, "EVALUATION_TAU");
  tau_gp_ = (tauloc == Inpar::ScaTra::evaltau_integration_point);  // set true/false

  const Inpar::ScaTra::EvalMat matloc =
      Core::UTILS::IntegralValue<Inpar::ScaTra::EvalMat>(stablist, "EVALUATION_MAT");
  mat_gp_ = (matloc == Inpar::ScaTra::evalmat_integration_point);  // set true/false

  // check for illegal combinations
  if (sgvel_ or assgd_)
  {
    // check for matching flags
    if (not mat_gp_ or not tau_gp_)
      FOUR_C_THROW(
          "Evaluation of material and stabilization parameters need to be done at the integration "
          "points if subgrid-scale velocity is included!");
  }

  // get quantities for finite difference check
  fdcheck_ = Core::UTILS::GetAsEnum<Inpar::ScaTra::FdCheck>(parameters, "fdcheck");
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

  // set flag for external force
  has_external_force_ = parameters.get<bool>("has_external_force", false);
}

int Discret::ELEMENTS::ScaTraEleParameterStd::nds_disp() const
{
  FOUR_C_ASSERT(nds_disp_ != -1,
      "You try to access the number of dofset associated with displacement dofs without "
      "having set it!");
  return nds_disp_;
}

int Discret::ELEMENTS::ScaTraEleParameterStd::nds_growth() const
{
  FOUR_C_ASSERT(nds_growth_ != -1,
      "You try to access the number of dofset associated with interface growth dofs without "
      "having set it!");
  return nds_growth_;
}

int Discret::ELEMENTS::ScaTraEleParameterStd::nds_micro() const
{
  FOUR_C_ASSERT(nds_micro_ != -1,
      "You try to access the number of dofset to write micro scale values on without having "
      "set it!");
  return nds_micro_;
}

int Discret::ELEMENTS::ScaTraEleParameterStd::nds_pres() const
{
  FOUR_C_ASSERT(nds_pres_ != -1,
      "You try to access the number of dofset associated with pressure dofs without having "
      "set it!");
  return nds_pres_;
}

int Discret::ELEMENTS::ScaTraEleParameterStd::nds_sca_tra() const
{
  FOUR_C_ASSERT(nds_scatra_ != -1,
      "You try to access the number of dofset associated with scalar transport dofs without having "
      "set it!");
  return nds_scatra_;
}

int Discret::ELEMENTS::ScaTraEleParameterStd::nds_thermo() const
{
  FOUR_C_ASSERT(nds_thermo_ != -1,
      "You try to access the number of dofset associated with temperature dofs without having "
      "set it!");
  return nds_thermo_;
}

int Discret::ELEMENTS::ScaTraEleParameterStd::nds_two_tensor_quantity() const
{
  FOUR_C_ASSERT(nds_two_tensor_quantitiy_ != -1,
      "You try to access the number of dofset associated with two-tensor quantity dofs without "
      "having set it!");
  return nds_two_tensor_quantitiy_;
}

int Discret::ELEMENTS::ScaTraEleParameterStd::nds_vel() const
{
  FOUR_C_ASSERT(nds_vel_ != -1,
      "You try to access the number of dofset associated with velocity related dofs without "
      "having set it!");
  return nds_vel_;
}

int Discret::ELEMENTS::ScaTraEleParameterStd::nds_wss() const
{
  FOUR_C_ASSERT(nds_wss_ != -1,
      "You try to access the number of dofset associated with wall shear stress dofs without "
      "having set it!");
  return nds_wss_;
}
FOUR_C_NAMESPACE_CLOSE
