/*----------------------------------------------------------------------*/
/*!
\file inpar_boolify.H

\brief Auxiliar routine to boolify integral Yes/No data

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifdef CCADISCRET

#ifndef INPAR_BOOLIFY_H
#define INPAR_BOOLIFY_H

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*/
/* content */
namespace DRT {

  namespace INPUT {

    /// Auxiliar routine to boolify integral Yes/No data
    ///
    /// Parameters consisting integral values for Yes/No tuples
    /// are removed and replaced by a bool parameter holding true/false
    /// respectively.
    void BoolifyValidInputParameters(
      Teuchos::ParameterList& list  ///< the valid input parameter list
      );

  }  // namespace INPUT
}  // namespace NOX

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_BOOLIFY_H
#endif  // #ifdef CCADISCRET
