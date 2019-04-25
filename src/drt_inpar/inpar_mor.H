/*----------------------------------------------------------------------*/
/*!
\file inpar_mor.H

\brief Input parameters for model order reduction

\level 2

\maintainer Amadeus Gebauer
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_MOR_H
#define INPAR_MOR_H

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
  namespace MOR
  {
    /// Defines all valid parameters for model order reduction
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace MOR
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MOR_H
