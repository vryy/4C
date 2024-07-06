/*---------------------------------------------------------------------*/
/*! \file
\brief singleton class holding all interface parameters required for boundary element evaluation

\level 2

*/
/*---------------------------------------------------------------------*/

#include "4C_scatra_ele_parameter_boundary.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraEleParameterBoundary*
Discret::ELEMENTS::ScaTraEleParameterBoundary::instance(const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const std::string& disname) {
        return std::unique_ptr<ScaTraEleParameterBoundary>(new ScaTraEleParameterBoundary(disname));
      });

  return singleton_map[disname].instance(Core::UTILS::SingletonAction::create, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraEleParameterBoundary::ScaTraEleParameterBoundary(
    const std::string& disname)
    : alphaa_(0.0),
      alphac_(0.0),
      conditiontype_(Core::Conditions::ConditionType::none),
      convtolimplicit_bv_(-1.0),
      density_(-1.0),
      molar_heat_capacity_(-1.0),
      is_pseudo_contact_(false),
      itemaxmimplicit_bv_(-1),
      kineticmodel_(-1),
      kr_(-1.0),
      molarmass_(-1.0),
      numelectrons_(0),
      numscal_(-1),
      onoff_(nullptr),
      peltier_(0.0),
      permeabilities_(nullptr),
      regularizationparameter_(-1.0),
      regularizationtype_(Inpar::S2I::RegularizationType::regularization_undefined),
      resistance_(0.0),
      resistivity_(0.0),
      capacitance_(0.0),
      stoichiometries_(nullptr),
      thermoperm_(-1.0)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_parameters(
    Teuchos::ParameterList& parameters)
{
  kineticmodel_ = parameters.get<int>("kinetic model", std::numeric_limits<int>::infinity());
  conditiontype_ = parameters.get<Core::Conditions::ConditionType>(
      "condition type", Core::Conditions::ConditionType::none);

  // set parameters to internal members depending on condition type
  switch (conditiontype_)
  {
    case Core::Conditions::ConditionType::S2IKinetics:
    {
      // set parameters to internal members depending on kinetic model
      switch (kineticmodel_)
      {
        case Inpar::S2I::kinetics_constperm:
        case Inpar::S2I::kinetics_linearperm:
        {
          set_is_pseudo_contact(parameters);
          set_num_scal(parameters);
          set_permeabilities(parameters);
          break;
        }

        case Inpar::S2I::kinetics_constantinterfaceresistance:
        {
          set_is_pseudo_contact(parameters);
          set_resistance(parameters);
          set_num_electrons(parameters);
          set_on_off(parameters);
          break;
        }

        case Inpar::S2I::kinetics_nointerfaceflux:
        {
          // do nothing
          break;
        }

        case Inpar::S2I::kinetics_butlervolmer:
        case Inpar::S2I::kinetics_butlervolmerlinearized:
        case Inpar::S2I::kinetics_butlervolmerreduced:
        case Inpar::S2I::kinetics_butlervolmerreducedcapacitance:
        case Inpar::S2I::kinetics_butlervolmerreducedlinearized:
        case Inpar::S2I::kinetics_butlervolmerpeltier:
        case Inpar::S2I::kinetics_butlervolmerresistance:
        case Inpar::S2I::kinetics_butlervolmerreducedthermoresistance:
        case Inpar::S2I::kinetics_butlervolmerreducedresistance:
        {
          set_alpha(parameters);
          set_charge_transfer_constant(parameters);
          set_is_pseudo_contact(parameters);
          set_num_electrons(parameters);
          set_num_scal(parameters);
          set_stoichiometries(parameters);
          if (kineticmodel_ == Inpar::S2I::kinetics_butlervolmerreducedcapacitance)
          {
            set_capacitance(parameters);
          }
          if (kineticmodel_ == Inpar::S2I::kinetics_butlervolmerpeltier)
            set_peltier(parameters);
          else if (kineticmodel_ == Inpar::S2I::kinetics_butlervolmerresistance or
                   kineticmodel_ == Inpar::S2I::kinetics_butlervolmerreducedresistance)
          {
            set_conv_tol_iter_num(parameters);
            set_resistance(parameters);
          }
          if (kineticmodel_ == Inpar::S2I::kinetics_butlervolmerreducedthermoresistance)
          {
            set_energy_substance_ratio(parameters);
            set_thermo_perm(parameters);
          }
          break;
        }

        default:
        {
          FOUR_C_THROW("Not implemented for this kinetic model: %i", kineticmodel_);
          break;
        }
      }

      // regularization is not relevant for scatra-scatra interface coupling without growth
      regularizationtype_ = Inpar::S2I::RegularizationType::regularization_none;

      break;
    }

    case Core::Conditions::ConditionType::S2IKineticsGrowth:
    {
      // set parameters to internal members depending on kinetic model
      switch (kineticmodel_)
      {
        case Inpar::S2I::growth_kinetics_butlervolmer:
        {
          set_alpha(parameters);
          set_charge_transfer_constant(parameters);
          set_density_molar_mass(parameters);
          set_num_electrons(parameters);
          set_num_scal(parameters);
          set_regularization(parameters);
          set_resistivity(parameters);
          set_stoichiometries(parameters);
          break;
        }

        default:
        {
          FOUR_C_THROW("Not implemented for this kinetic model: %i", kineticmodel_);
          break;
        }
      }
      break;
    }

    default:
    {
      FOUR_C_THROW("Not implemented for this condition type: %i", conditiontype_);
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_alpha(Teuchos::ParameterList& parameters)
{
  alphaa_ = parameters.get<double>("alpha_a", std::numeric_limits<double>::infinity());
  alphac_ = parameters.get<double>("alpha_c", std::numeric_limits<double>::infinity());
  if (alphaa_ <= 0.0) FOUR_C_THROW("Alpha a must be positive!");
  if (alphac_ <= 0.0) FOUR_C_THROW("Alpha c must be positive!");
  if (alphaa_ + alphac_ != 1.0) FOUR_C_THROW("Sum of Alpha a and Alpha c must be equal to one!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_charge_transfer_constant(
    Teuchos::ParameterList& parameters)
{
  kr_ = parameters.get<double>("k_r", -1.0);
  if (kr_ <= 0.0) FOUR_C_THROW("Charge transfer constant k_r is negative!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_conv_tol_iter_num(
    Teuchos::ParameterList& parameters)
{
  convtolimplicit_bv_ =
      parameters.get<double>("CONVTOL_IMPLBUTLERVOLMER", std::numeric_limits<double>::infinity());
  if (convtolimplicit_bv_ <= 0.0) FOUR_C_THROW("Tolerance of convergence must be positive!");
  itemaxmimplicit_bv_ =
      parameters.get<int>("ITEMAX_IMPLBUTLERVOLMER", std::numeric_limits<int>::infinity());
  if (itemaxmimplicit_bv_ <= 0) FOUR_C_THROW("Maximum number of iterations must be positive!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_density_molar_mass(
    Teuchos::ParameterList& parameters)
{
  density_ = parameters.get<double>("density", std::numeric_limits<double>::infinity());
  if (density_ <= 0.0) FOUR_C_THROW("Density must be positive");

  molarmass_ = parameters.get<double>("molar mass", std::numeric_limits<double>::infinity());
  if (molarmass_ <= 0.0) FOUR_C_THROW("Molar mass must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_energy_substance_ratio(
    Teuchos::ParameterList& parameters)
{
  molar_heat_capacity_ =
      parameters.get<double>("molar_heat_capacity", std::numeric_limits<double>::infinity());
  if (molar_heat_capacity_ < 0.0) FOUR_C_THROW("Ratio of energy- and mass-flux must be positive!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_is_pseudo_contact(
    Teuchos::ParameterList& parameters)
{
  is_pseudo_contact_ =
      (parameters.get<int>("is_pseudo_contact", std::numeric_limits<int>::infinity()) == 1);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_num_electrons(
    Teuchos::ParameterList& parameters)
{
  numelectrons_ = parameters.get<int>("numelectrons", std::numeric_limits<int>::infinity());
  if (numelectrons_ != 1)
    FOUR_C_THROW(
        "Invalid number of electrons in charge transfer at electrode-electrolyte interface!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_num_scal(Teuchos::ParameterList& parameters)
{
  numscal_ = parameters.get<int>("numscal", std::numeric_limits<int>::infinity());
  if (numscal_ <= 0) FOUR_C_THROW("Scalar must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_peltier(Teuchos::ParameterList& parameters)
{
  peltier_ = parameters.get<double>("peltier", std::numeric_limits<double>::infinity());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_permeabilities(
    Teuchos::ParameterList& parameters)
{
  permeabilities_ = parameters.get<const std::vector<double>*>("permeabilities");
  for (auto permeability : *permeabilities_)
    if (permeability < 0.0) FOUR_C_THROW("Permeability must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_regularization(
    Teuchos::ParameterList& parameters)
{
  regularizationparameter_ = parameters.get<double>("regpar", -1.0);
  if (regularizationparameter_ < 0.0)
    FOUR_C_THROW("Regularization parameter for lithium stripping must not be negative!");
  regularizationtype_ = static_cast<Inpar::S2I::RegularizationType>(
      parameters.get<int>("regtype", std::numeric_limits<int>::infinity()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_resistance(
    Teuchos::ParameterList& parameters)
{
  resistance_ = parameters.get<double>("resistance", std::numeric_limits<double>::infinity());
  if (resistance_ <= 0.0) FOUR_C_THROW("Resistance must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_resistivity(
    Teuchos::ParameterList& parameters)
{
  resistivity_ = 1.0 / (parameters.get<double>("conductivity", -1.0));
  if (resistivity_ <= 0.0) FOUR_C_THROW("Conductivity must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_capacitance(
    Teuchos::ParameterList& parameters)
{
  capacitance_ = parameters.get<double>("capacitance", -1.0);
  if (capacitance_ <= 0.0) FOUR_C_THROW("Capacitance must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_stoichiometries(
    Teuchos::ParameterList& parameters)
{
  stoichiometries_ = parameters.get<const std::vector<int>*>("stoichiometries");

  if (stoichiometries_ == nullptr)
  {
    FOUR_C_THROW(
        "Cannot get vector of stoichiometric coefficients for scatra-scatra interface coupling!");
  }
  else
  {
    if (stoichiometries_->size() != 1)
      FOUR_C_THROW("Number of stoichiometric coefficients does not match number of scalars!");
    if ((*stoichiometries_)[0] != -1) FOUR_C_THROW("Invalid stoichiometric coefficient!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_thermo_perm(
    Teuchos::ParameterList& parameters)
{
  thermoperm_ = parameters.get<double>("thermoperm", std::numeric_limits<double>::infinity());
  if (thermoperm_ <= 0.0) FOUR_C_THROW("Thermo permeability must be posititve!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterBoundary::set_on_off(Teuchos::ParameterList& parameters)
{
  onoff_ = parameters.get<const std::vector<int>*>("onoff");
  if (onoff_ == nullptr) FOUR_C_THROW("Cannot get vector 'onoff' from parameter list");
  if (onoff_->size() != 2) FOUR_C_THROW("Only two dofs are supported");
}

FOUR_C_NAMESPACE_CLOSE
