/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for model order reduction

\level 2

*/

/*----------------------------------------------------------------------*/
#ifndef BACI_INPAR_MOR_H
#define BACI_INPAR_MOR_H

#include "baci_inpar_parameterlist_utils.H"

BACI_NAMESPACE_OPEN

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
BACI_NAMESPACE_CLOSE

#endif  // INPAR_MOR_H
