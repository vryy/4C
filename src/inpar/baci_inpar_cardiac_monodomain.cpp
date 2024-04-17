/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for cardiac monodomain

\level 2

*/
/*----------------------------------------------------------------------*/
#include "baci_inpar_cardiac_monodomain.hpp"

#include "baci_lib_conditiondefinition.hpp"
#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void INPAR::EP::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& epcontrol = list->sublist("CARDIAC MONODOMAIN CONTROL", false,
      "control parameters for cardiac electrophysiology problems\n");


  // Parameters for reaction-diffusion systems (for example cardiac electrophysiology)
  CORE::UTILS::IntParameter("WRITEMAXINTSTATE", 0,
      "number of maximal internal state variables to be postprocessed", &epcontrol);
  CORE::UTILS::IntParameter("WRITEMAXIONICCURRENTS", 0,
      "number of maximal ionic currents to be postprocessed", &epcontrol);

  CORE::UTILS::DoubleParameter("ACTTHRES", 1.0,
      "threshold for the potential for computing and postprocessing activation time ", &epcontrol);
}


void INPAR::EP::SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist)
{
}

FOUR_C_NAMESPACE_CLOSE
