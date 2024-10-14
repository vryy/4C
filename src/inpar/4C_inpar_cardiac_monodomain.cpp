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

void Inpar::ElectroPhysiology::set_valid_parameters(Teuchos::ParameterList& list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& epcontrol = list.sublist("CARDIAC MONODOMAIN CONTROL", false,
      "control parameters for cardiac electrophysiology problems\n");


  // Parameters for reaction-diffusion systems (for example cardiac electrophysiology)
  Core::Utils::int_parameter("WRITEMAXINTSTATE", 0,
      "number of maximal internal state variables to be postprocessed", &epcontrol);
  Core::Utils::int_parameter("WRITEMAXIONICCURRENTS", 0,
      "number of maximal ionic currents to be postprocessed", &epcontrol);

  Core::Utils::double_parameter("ACTTHRES", 1.0,
      "threshold for the potential for computing and postprocessing activation time ", &epcontrol);
}


void Inpar::ElectroPhysiology::set_valid_conditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
}

FOUR_C_NAMESPACE_CLOSE
