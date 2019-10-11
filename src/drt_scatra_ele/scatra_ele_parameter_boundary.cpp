/*---------------------------------------------------------------------*/
/*! \file
\brief singleton class holding all interface parameters required for boundary element evaluation

\level 2

\maintainer Christoph Schmidt
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
  if (delete_me == NULL)
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
        return NULL;
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

  return;
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                   civaner 08/19 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterBoundary::ScaTraEleParameterBoundary(const std::string& disname)
    : alphaa_(0.0),
      alphac_(0.0),
      conditiontype_(DRT::Condition::ConditionType::none),
      density_(-1.0),
      kineticmodel_(-1),
      kr_(-1.0),
      molarmass_(-1.0),
      numelectrons_(0),
      numscal_(-1),
      peltier_(0.0),
      permeabilities_(NULL),
      regularizationparameter_(-1.0),
      regularizationtype_("undefined"),
      resistivity_(0.0),
      stoichiometries_(NULL),
      resistance_(0.0),
      convtolmodifiedBV_(-1.0),
      itemaxmodifiedBV_(-1.0)
{
  return;
}

/*----------------------------------------------------------------------*
 | set parameters                                         civaner 08/19 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterBoundary::SetParameters(Teuchos::ParameterList& parameters)
{
  kineticmodel_ = parameters.get<int>("kinetic model", std::numeric_limits<int>::infinity());
  numscal_ = parameters.get<int>("numscal", std::numeric_limits<int>::infinity());
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
          permeabilities_ = parameters.get<std::vector<double>*>("permeabilities");
          break;
        }

        case INPAR::S2I::kinetics_butlervolmer:
        case INPAR::S2I::kinetics_butlervolmerreduced:
        case INPAR::S2I::kinetics_butlervolmerpeltier:
        case INPAR::S2I::kinetics_butlervolmerresistance:
        {
          stoichiometries_ = parameters.get<std::vector<int>*>("stoichiometries");
          numelectrons_ = parameters.get<int>("numelectrons", std::numeric_limits<int>::infinity());
          kr_ = parameters.get<double>("k_r", -1.0);
          alphaa_ = parameters.get<double>("alpha_a", std::numeric_limits<double>::infinity());
          alphac_ = parameters.get<double>("alpha_c", std::numeric_limits<double>::infinity());

          if (kineticmodel_ == INPAR::S2I::kinetics_butlervolmerpeltier)
            peltier_ = parameters.get<double>("peltier", std::numeric_limits<double>::infinity());

          if (kineticmodel_ == INPAR::S2I::kinetics_butlervolmerresistance)
          {
            resistance_ =
                parameters.get<double>("resistance", std::numeric_limits<double>::infinity());
            convtolmodifiedBV_ = parameters.get<double>(
                "CONVTOL_MODBUTLERVOLMER", std::numeric_limits<double>::infinity());
            itemaxmodifiedBV_ = parameters.get<double>(
                "ITEMAX_MODBUTLERVOLMER", std::numeric_limits<double>::infinity());
          }
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

    case DRT::Condition::ConditionType::S2ICouplingGrowth:
    {
      // set parameters to internal members depending on kinetic model
      switch (kineticmodel_)
      {
        case INPAR::S2I::growth_kinetics_butlervolmer:
        {
          stoichiometries_ = parameters.get<std::vector<int>*>("stoichiometries");
          numelectrons_ = parameters.get<int>("numelectrons", std::numeric_limits<int>::infinity());
          kr_ = parameters.get<double>("k_r", -1.0);
          alphaa_ = parameters.get<double>("alpha_a", std::numeric_limits<double>::infinity());
          alphac_ = parameters.get<double>("alpha_c", std::numeric_limits<double>::infinity());
          density_ = parameters.get<double>("density", std::numeric_limits<double>::infinity());
          molarmass_ =
              parameters.get<double>("molar mass", std::numeric_limits<double>::infinity());
          regularizationparameter_ = parameters.get<double>("regpar", -1.0);
          regularizationtype_ = parameters.get<std::string>("regtype");
          resistivity_ = 1.0 / (parameters.get<double>("conductivity", -1.0));

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

  return;
}
