/*----------------------------------------------------------------------*/
/*! \file


\brief Input parameters for model order reduction

\level 2
*/

/*----------------------------------------------------------------------*/

#include "baci_inpar_mor.H"

#include "baci_inpar_validparameters.H"

BACI_NAMESPACE_OPEN

void INPAR::MOR::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& mor = list->sublist("MOR", false, "");

  StringParameter("POD_MATRIX", "none", "filename of file containing projection matrix", &mor);
}

BACI_NAMESPACE_CLOSE
