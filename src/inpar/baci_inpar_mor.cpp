/*----------------------------------------------------------------------*/
/*! \file


\brief Input parameters for model order reduction

\level 2
*/

/*----------------------------------------------------------------------*/

#include "baci_inpar_mor.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

void INPAR::MOR::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& mor = list->sublist("MOR", false, "");

  CORE::UTILS::StringParameter(
      "POD_MATRIX", "none", "filename of file containing projection matrix", &mor);
}

BACI_NAMESPACE_CLOSE
