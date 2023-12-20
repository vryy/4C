/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-interaction

\level 2


 *------------------------------------------------------------------------------------------------*/


#ifndef BACI_INPAR_SSTI_H
#define BACI_INPAR_SSTI_H

#include "baci_config.H"

#include "baci_inpar_parameterlist_utils.H"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

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
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace SSTI
}  // namespace INPAR

BACI_NAMESPACE_CLOSE

#endif  // INPAR_SSTI_H
