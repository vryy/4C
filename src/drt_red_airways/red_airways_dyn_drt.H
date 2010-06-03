/*!----------------------------------------------------------------------
\file red_airways_dyn_drt.H
\brief Control routine for reduced dimesional airways network
 (time) integration. Calls


<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
/************/
#ifdef CCADISCRET
#ifndef RED_AIRWAYS_DRT_DRT_H
#define RED_AIRWAYS_DRT_DRT_H

#ifdef __cplusplus


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_io/io.H"

using namespace std;
using namespace Teuchos;

/* headers and prototypes of functions also accessible from c */
extern "C"
{
  void dyn_red_airways_drt();
}

#endif /* endif __cplusplus*/

#endif  /* #ifndef RED_AIRWAYS_DRT_DRT_H */
#endif  /* #ifdef CCADISCRET */

