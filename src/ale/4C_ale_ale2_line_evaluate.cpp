/*----------------------------------------------------------------------------*/
/*! \file

\brief Evaluate ALE2 elements

\level 1

*/
/*----------------------------------------------------------------------------*/

#include "4C_ale_ale2.hpp"
#include "4C_lib_discret.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale2Line::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  FOUR_C_THROW("Line Neumann condition not implemented for Ale2");

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
