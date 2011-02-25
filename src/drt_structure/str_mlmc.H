/*----------------------------------------------------------------------*/
/*!
\file str_mlmc.H
\brief

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef STR_MLMC_H
#define STR_MLMC_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*/
/* Structural Dynamics */
namespace STR
{

  /*====================================================================*/
  //! Multi Level Monte Caro analysis of structures
  void mlmc();

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  /* #ifndef STRUDYN_DIRECT_H */
#endif  /* #ifdef CCADISCRET */
