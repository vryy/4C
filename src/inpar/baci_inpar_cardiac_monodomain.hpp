/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for cardiac monodomain

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_INPAR_CARDIAC_MONODOMAIN_HPP
#define BACI_INPAR_CARDIAC_MONODOMAIN_HPP


#include "baci_config.hpp"

#include "baci_inpar_parameterlist_utils.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

namespace INPAR
{
  namespace EP
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
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);
  }  // namespace EP
}  // namespace INPAR
BACI_NAMESPACE_CLOSE

#endif  // INPAR_CARDIAC_MONODOMAIN_H
