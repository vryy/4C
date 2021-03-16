/*----------------------------------------------------------------------*/
/*! \file

\brief parameter class for elch on manifold kinetics

\level 2

*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_dserror.H"

#include "scatra_ele_parameter_elch_manifold.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElchManifold*
DRT::ELEMENTS::ScaTraEleParameterElchManifold::Instance(
    const std::string& disname, const ScaTraEleParameterElchManifold* delete_me)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, ScaTraEleParameterElchManifold*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == nullptr)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterElchManifold(disname);
  }

  // destruct instance
  else
  {
    for (auto instance = instances.begin(); instance != instances.end(); ++instance)
    {
      if (instance->second == delete_me)
      {
        delete instance->second;
        instances.erase(instance);
        return nullptr;
      }
    }
    dserror("Could not locate the desired instance. Internal error.");
  }

  // return existing or newly created instance
  return instances[disname];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElchManifold::Done()
{
  // delete singleton
  Instance("", this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElchManifold::ScaTraEleParameterElchManifold(
    const std::string& disname)
    : kinetic_model_(INPAR::SSI::kinetics_noflux),
      num_electrons_(-1),
      resistance_(-1.0),
      use_other_side_(false)

{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElchManifold::SetParameters(
    Teuchos::ParameterList& parameters)
{
  kinetic_model_ = parameters.get<int>("kinetic_model");

  if (kinetic_model_ == INPAR::SSI::kinetics_constantinterfaceresistance)
  {
    resistance_ = parameters.get<double>("resistance");

    num_electrons_ = parameters.get<int>("num_electrons");
  }

  use_other_side_ = parameters.get<bool>("use_other_side");

  // safety checks
  if (kinetic_model_ == INPAR::SSI::kinetics_constantinterfaceresistance)
  {
    if (resistance_ <= 0.0) dserror("Resistance is non-positive!");
    if (num_electrons_ <= 0) dserror("Number of electrons must be positive!");
  }
}
