/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for model order reduction

\level 2

*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_MOR_H
#define INPAR_MOR_H

#include "inpar_parameterlist_utils.H"

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
