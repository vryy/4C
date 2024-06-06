/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-interaction

\level 2


 *------------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_INPAR_SSTI_HPP
#define FOUR_C_INPAR_SSTI_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
namespace Inpar
{
  namespace SSTI
  {
    /// Type of coupling strategy for SSI problems
    enum class SolutionScheme
    {
      monolithic
    };

    //! type of scalar transport time integration
    enum class ScaTraTimIntType
    {
      elch
    };

    /// set the ssti parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific ssti conditions
    void SetValidConditions(
        std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist);

  }  // namespace SSTI
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

#endif
