/*----------------------------------------------------------------------------*/
/*! \file

\brief Evaluate ALE2 elements

\level 1

*/
/*----------------------------------------------------------------------------*/

#include "4C_ale_ale2.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int Discret::ELEMENTS::Ale2Line::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  FOUR_C_THROW("Line Neumann condition not implemented for Ale2");

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
