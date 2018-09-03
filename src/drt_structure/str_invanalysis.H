/*----------------------------------------------------------------------*/
/*!
\file str_invanalysis.H
\brief Inverse analysis for structures

\level 1

<pre>
\maintainer Harald Willmann
            willmann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef STR_INVANALYSIS_H
#define STR_INVANALYSIS_H

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
  //! Inverse analysis of structures
  //!
  //! \author rausch \date 02/09
  void invanalysis();

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif /* #ifndef STRUDYN_DIRECT_H */
