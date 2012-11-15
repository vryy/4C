/*----------------------------------------------------------------------*/
/*!
\file str_mlmc.H
\brief

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef STR_MLMC_H
#define STR_MLMC_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Time.h"
#include "Teuchos_RCP.hpp"


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
