/*----------------------------------------------------------------------------*/
/*! \file

\brief Evaluate ALE2 elements

\level 1

*/
/*----------------------------------------------------------------------------*/

#include "baci_ale_ale2.H"
#include "baci_lib_discret.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale2Line::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  dserror("Line Neumann condition not implemented for Ale2");

  return 0;
}

BACI_NAMESPACE_CLOSE
