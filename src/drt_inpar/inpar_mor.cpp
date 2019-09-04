/*----------------------------------------------------------------------*/
/*! \file

\maintainer Amadeus Gebauer

\brief Input parameters for model order reduction

\level 2
*/

/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_mor.H"

void INPAR::MOR::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& mor = list->sublist("MOR", false, "");

  StringParameter("POD_MATRIX", "none", "filename of file containing projection matrix", &mor);
}
