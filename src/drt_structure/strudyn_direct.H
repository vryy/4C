/*----------------------------------------------------------------------*/
/*!
\file strudyn_direct.H
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef STRUDYN_DIRECT_H
#define STRUDYN_DIRECT_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"

/*----------------------------------------------------------------------*/
/* Structural Dynamics */
namespace STR
{

  /*====================================================================*/
  //! Direct time integration of temporally (and spatially) discretised
  //! structure
  //!
  //! \author bborn \date 06/08
  void strudyn_direct();

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  /* #ifndef STRUDYN_DIRECT_H */
#endif  /* #ifdef CCADISCRET */
