/*----------------------------------------------------------------------*/
/*! \file

\brief Auxiliar routine to boolify integral Yes/No data


\level 0
*/

/*----------------------------------------------------------------------*/
/* definitions */

#ifndef INPAR_BOOLIFYPARAMETERS_H
#define INPAR_BOOLIFYPARAMETERS_H

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*/
/* content */
namespace DRT
{
  namespace INPUT
  {
    /// Auxiliar routine to boolify integral Yes/No data
    ///
    /// Parameters consisting integral values for Yes/No tuples
    /// are removed and replaced by a bool parameter holding true/false
    /// respectively.
    void BoolifyValidInputParameters(
        Teuchos::ParameterList& list  ///< the valid input parameter list
    );

  }  // namespace INPUT
}  // namespace DRT

/*----------------------------------------------------------------------*/
#endif  // INPAR_BOOLIFYPARAMETERS_H
