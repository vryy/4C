/*----------------------------------------------------------------------*/
/*!
\file drt_boolifyparameters.H

\brief Auxiliar routine to boolify integral Yes/No data

\maintainer Martin Kronbichler

\level 0
*/

/*----------------------------------------------------------------------*/
/* definitions */

#ifndef DRT_BOOLIFYPARAMETERS_H
#define DRT_BOOLIFYPARAMETERS_H

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
#endif  // #ifndef DRT_BOOLIFYPARAMETERS_H
