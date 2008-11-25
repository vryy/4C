/*----------------------------------------------------------------------*/
/*!
\file str_aux.H
\brief Various auxiliar methods needed in structural analysis

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifdef CCADISCRET
#ifndef STR_UTILS_H
#define STR_UTILS_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"

#include "../drt_inpar/inpar_structure.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /// A namespace for the various unsorted auxiliar kind of
  /// methods needed in structural analysis
  ///
  /// \author bborn
  /// \date 11/08
  namespace AUX
  {

    /// Determine norm of force residual
    const double CalculateVectorNorm
    (
      const enum INPAR::STR::VectorNorm norm,  ///< norm to use
      const Teuchos::RCP<Epetra_Vector> vect  ///< the vector of interest
    );

  }  // namespace AUX

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STR_UTILS_H
#endif  // #ifdef CCADISCRET
