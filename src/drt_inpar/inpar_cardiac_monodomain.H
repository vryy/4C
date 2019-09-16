/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for cardiac monodomain

\level 2

\maintainer Anh-Tu Vuong
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_CARDIAC_MONODOMAIN_H
#define INPAR_CARDIAC_MONODOMAIN_H


#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

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
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace EP
}  // namespace INPAR
#endif  // #ifndef INPAR_CARDIAC_MONODOMAIN_H
