/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for inverse analysis

\level 2


*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_INVANALYSIS_H
#define INPAR_INVANALYSIS_H

#include "inpar_parameterlist_utils.H"

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
