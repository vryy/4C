/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for cardiac monodomain

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_CARDIAC_MONODOMAIN_HPP
#define FOUR_C_INPAR_CARDIAC_MONODOMAIN_HPP


#include "4C_config.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration

namespace Inpar
{
  namespace ElectroPhysiology
  {
    /// possible types of evaluation of reaction term
    enum EvalType
    {
      ep_implicit,
      ep_semi_implicit,
    };

    /// set the elch parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific elch conditions
    void SetValidConditions(
        std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist);
  }  // namespace ElectroPhysiology
}  // namespace Inpar
FOUR_C_NAMESPACE_CLOSE

#endif
