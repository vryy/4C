/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for inverse analysis

\level 2


*/

/*----------------------------------------------------------------------*/
#ifndef BACI_INPAR_INVANALYSIS_H
#define BACI_INPAR_INVANALYSIS_H

#include "baci_inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace INVANA
  {
    /// set the invana parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific inverse analysis conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace INVANA

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // INPAR_INVANALYSIS_H
