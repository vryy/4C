/*---------------------------------------------------------------------*/
/*! \file
\brief singleton class holding all interface parameters required for boundary element evaluation

\level 2

*/
/*---------------------------------------------------------------------*/

#include "../drt_lib/drt_dserror.H"
#include "scatra_ele_parameter_boundary.H"
#include "../drt_inpar/inpar_s2i.H"

/*----------------------------------------------------------------------*
 | singleton access method                                civaner 08/19 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterBoundary* DRT::ELEMENTS::ScaTraEleParameterBoundary::Instance(
    const std::string& disname, const ScaTraEleParameterBoundary* delete_me)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, ScaTraEleParameterBoundary*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == nullptr)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterBoundary(disname);
  }

  // destruct instance
  else
  {
    for (std::map<std::string, ScaTraEleParameterBoundary*>::iterator i = instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return nullptr;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  // return existing or newly created instance
  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                  civaner 08/19 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::Done()
{
  // delete singleton
  Instance("", this);
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                   civaner 08/19 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterBoundary::ScaTraEleParameterBoundary(const std::string& disname)
    : alphaa_(0.0),
      alphac_(0.0),
      conditiontype_(DRT::Condition::ConditionType::none),
      convtolimplicitBV_(-1.0),
      density_(-1.0),
      itemaxmimplicitBV_(-1),
      kineticmodel_(-1),
      kr_(-1.0),
      molarmass_(-1.0),
      numelectrons_(0),
      numscal_(-1),
      peltier_(0.0),
      permeabilities_(nullptr),
      regularizationparameter_(-1.0),
      regularizationtype_(INPAR::S2I::RegularizationType::regularization_undefined),
      resistance_(0.0),
      resistivity_(0.0),
      stoichiometries_(nullptr)
{
  return;
}

/*----------------------------------------------------------------------*
 | set parameters                                         civaner 08/19 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetParameters(Teuchos::ParameterList& parameters)
{
  kineticmodel_ = parameters.get<int>("kinetic model", std::numeric_limits<int>::infinity());
  conditiontype_ = parameters.get<DRT::Condition::ConditionType>(
      "condition type", DRT::Condition::ConditionType::none);

  // set parameters to internal members depending on condition type
  switch (conditiontype_)
  {
    case DRT::Condition::ConditionType::S2ICoupling:
    {
      // set parameters to internal members depending on kinetic model
      switch (kineticmodel_)
      {
        case INPAR::S2I::kinetics_constperm:
        {
          SetNumScal(parameters);
          SetPermeabilities(parameters);
          break;
        }

        case INPAR::S2I::kinetics_constantinterfaceresistance:
        {
          SetResistance(parameters);
          break;
        }

        case INPAR::S2I::kinetics_nointerfaceflux:
        {
          // do nothing
          break;
        }

        case INPAR::S2I::kinetics_butlervolmer:
        case INPAR::S2I::kinetics_butlervolmerreduced:
        case INPAR::S2I::kinetics_butlervolmerpeltier:
        case INPAR::S2I::kinetics_butlervolmerresistance:
        case INPAR::S2I::kinetics_butlervolmerreducedwithresistance:
        {
          SetAlpha(parameters);
          SetChargeTransferConstant(parameters);
          SetNumElectrons(parameters);
          SetNumScal(parameters);
          SetStoichiometries(parameters);
          if (kineticmodel_ == INPAR::S2I::kinetics_butlervolmerpeltier)
            SetPeltier(parameters);
          else if (kineticmodel_ == INPAR::S2I::kinetics_butlervolmerresistance or
                   kineticmodel_ == INPAR::S2I::kinetics_butlervolmerreducedwithresistance)
          {
            SetConvTolIterNum(parameters);
            SetResistance(parameters);
          }
          break;
        }

        default:
        {
          dserror("Not implemented for this kinetic model: %i", kineticmodel_);
          break;
        }
      }

      // regularization is not relevant for scatra-scatra interface coupling without growth
      regularizationtype_ = INPAR::S2I::RegularizationType::regularization_none;

      break;
    }

    case DRT::Condition::ConditionType::S2ICouplingGrowth:
    {
      // set parameters to internal members depending on kinetic model
      switch (kineticmodel_)
      {
        case INPAR::S2I::growth_kinetics_butlervolmer:
        {
          SetAlpha(parameters);
          SetChargeTransferConstant(parameters);
          SetDensityMolarMass(parameters);
          SetNumElectrons(parameters);
          SetNumScal(parameters);
          SetRegularization(parameters);
          SetResistivity(parameters);
          SetStoichiometries(parameters);
          break;
        }

        default:
        {
          dserror("Not implemented for this kinetic model: %i", kineticmodel_);
          break;
        }
      }
      break;
    }

    default:
    {
      dserror("Not implemented for this condition type: %i", conditiontype_);
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetAlpha(Teuchos::ParameterList& parameters)
{
  alphaa_ = parameters.get<double>("alpha_a", std::numeric_limits<double>::infinity());
  alphac_ = parameters.get<double>("alpha_c", std::numeric_limits<double>::infinity());
  if (alphaa_ <= 0.0) dserror("Alpha a must be positive!");
  if (alphac_ <= 0.0) dserror("Alpha c must be positive!");
  if (alphaa_ + alphac_ != 1.0) dserror("Sum of Alpha a and Alpha c must be equal to one!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetChargeTransferConstant(
    Teuchos::ParameterList& parameters)
{
  kr_ = parameters.get<double>("k_r", -1.0);
  if (kr_ <= 0.0) dserror("Charge transfer constant k_r is negative!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetConvTolIterNum(
    Teuchos::ParameterList& parameters)
{
  convtolimplicitBV_ =
      parameters.get<double>("CONVTOL_IMPLBUTLERVOLMER", std::numeric_limits<double>::infinity());
  if (convtolimplicitBV_ <= 0.0) dserror("Tolerance of convergence must be positive!");
  itemaxmimplicitBV_ =
      parameters.get<int>("ITEMAX_IMPLBUTLERVOLMER", std::numeric_limits<int>::infinity());
  if (itemaxmimplicitBV_ <= 0) dserror("Maximum number of iterations must be positive!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetDensityMolarMass(
    Teuchos::ParameterList& parameters)
{
  density_ = parameters.get<double>("density", std::numeric_limits<double>::infinity());
  if (density_ <= 0.0) dserror("Density must be positive");

  molarmass_ = parameters.get<double>("molar mass", std::numeric_limits<double>::infinity());
  if (molarmass_ <= 0.0) dserror("Molar mass must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetNumElectrons(Teuchos::ParameterList& parameters)
{
  numelectrons_ = parameters.get<int>("numelectrons", std::numeric_limits<int>::infinity());
  if (numelectrons_ != 1)
    dserror("Invalid number of electrons in charge transfer at electrode-electrolyte interface!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetNumScal(Teuchos::ParameterList& parameters)
{
  numscal_ = parameters.get<int>("numscal", std::numeric_limits<int>::infinity());
  if (numscal_ <= 0) dserror("Scalar must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetPeltier(Teuchos::ParameterList& parameters)
{
  peltier_ = parameters.get<double>("peltier", std::numeric_limits<double>::infinity());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetPermeabilities(
    Teuchos::ParameterList& parameters)
{
  permeabilities_ = parameters.get<std::vector<double>*>("permeabilities");
  for (auto permeability : *permeabilities_)
    if (permeability <= 0.0) dserror("Permeability must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetRegularization(
    Teuchos::ParameterList& parameters)
{
  regularizationparameter_ = parameters.get<double>("regpar", -1.0);
  if (regularizationparameter_ < 0.0)
    dserror("Regularization parameter for lithium stripping must not be negative!");
  regularizationtype_ = static_cast<INPAR::S2I::RegularizationType>(
      parameters.get<int>("regtype", std::numeric_limits<int>::infinity()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetResistance(Teuchos::ParameterList& parameters)
{
  resistance_ = parameters.get<double>("resistance", std::numeric_limits<double>::infinity());
  if (resistance_ <= 0.0) dserror("Resistance must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetResistivity(Teuchos::ParameterList& parameters)
{
  resistivity_ = 1.0 / (parameters.get<double>("conductivity", -1.0));
  if (resistivity_ <= 0.0) dserror("Conductivity must be positive");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetStoichiometries(
    Teuchos::ParameterList& parameters)
{
  stoichiometries_ = parameters.get<std::vector<int>*>("stoichiometries");

  if (stoichiometries_ == nullptr)
    dserror(
        "Cannot get vector of stoichiometric coefficients for scatra-scatra interface coupling!");
  if (stoichiometries_->size() != 1)
    dserror("Number of stoichiometric coefficients does not match number of scalars!");
  if ((*stoichiometries_)[0] != -1) dserror("Invalid stoichiometric coefficient!");
}
