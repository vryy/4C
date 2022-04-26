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
DRT::ELEMENTS::ScaTraEleParameterElchManifold::Instance(const std::string& disname)
{
  auto& selected = instances_[disname];

  // check whether instance already exists for current discretization, otherwise create it
  if (selected == nullptr)
    selected = std::unique_ptr<ScaTraEleParameterElchManifold>(
        new ScaTraEleParameterElchManifold(disname));

  return selected.get();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElchManifold::Done()
{
  for (const auto& instance : instances_)
  {
    if (instance.second.get() == this)
    {
      instances_.erase(instance.first);
      break;
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElchManifold::ScaTraEleParameterElchManifold(
    const std::string& disname)
    : evaluate_conc_flux_(false),
      evaluate_master_side_(false),
      evaluate_pot_flux_(false),
      kinetic_model_(INPAR::SSI::kinetics_noflux),
      num_electrons_(-1),
      resistance_(-1.0)
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
    if (resistance_ <= 0.0) dserror("Resistance is non-positive!");

    num_electrons_ = parameters.get<int>("num_electrons");
    if (num_electrons_ <= 0) dserror("Number of electrons must be positive!");

    evaluate_conc_flux_ = parameters.get<bool>("conc_flux");

    evaluate_pot_flux_ = parameters.get<bool>("pot_flux");
  }

  evaluate_master_side_ = parameters.get<bool>("evaluate_master_side");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<std::string, std::unique_ptr<DRT::ELEMENTS::ScaTraEleParameterElchManifold>>
    DRT::ELEMENTS::ScaTraEleParameterElchManifold::instances_;