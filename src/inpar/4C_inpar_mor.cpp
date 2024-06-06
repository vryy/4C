/*----------------------------------------------------------------------*/
/*! \file


\brief Input parameters for model order reduction

\level 2
*/

/*----------------------------------------------------------------------*/

#include "4C_inpar_mor.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void Inpar::ModelOrderRed::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& mor = list->sublist("MOR", false, "");

  Core::UTILS::StringParameter(
      "POD_MATRIX", "none", "filename of file containing projection matrix", &mor);
}

FOUR_C_NAMESPACE_CLOSE
