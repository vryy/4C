/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for cardiac monodomain

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_inpar_cardiac_monodomain.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void Inpar::ElectroPhysiology::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& epcontrol = list->sublist("CARDIAC MONODOMAIN CONTROL", false,
      "control parameters for cardiac electrophysiology problems\n");


  // Parameters for reaction-diffusion systems (for example cardiac electrophysiology)
  Core::UTILS::IntParameter("WRITEMAXINTSTATE", 0,
      "number of maximal internal state variables to be postprocessed", &epcontrol);
  Core::UTILS::IntParameter("WRITEMAXIONICCURRENTS", 0,
      "number of maximal ionic currents to be postprocessed", &epcontrol);

  Core::UTILS::DoubleParameter("ACTTHRES", 1.0,
      "threshold for the potential for computing and postprocessing activation time ", &epcontrol);
}


void Inpar::ElectroPhysiology::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
}

FOUR_C_NAMESPACE_CLOSE
