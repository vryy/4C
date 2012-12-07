/*----------------------------------------------------------------------*/
/*!
\file str_statinvanalysis.H
\brief

<pre>
Maintainer: Sebastian Kehl
            Kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de

</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef STR_STATINVANALYSIS_H
#define STR_STATINVANALYSIS_H

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
  //! Statistical inverse analysis of structures

  void statinvanalysis();

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  /* #ifndef STRUDYN_DIRECT_H */
