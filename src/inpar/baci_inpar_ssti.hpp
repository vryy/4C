/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-interaction

\level 2


 *------------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_INPAR_SSTI_HPP
#define FOUR_C_INPAR_SSTI_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

namespace INPUT
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
namespace INPAR
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
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace SSTI
}  // namespace INPAR

FOUR_C_NAMESPACE_CLOSE

#endif
