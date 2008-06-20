/*----------------------------------------------------------------------*/
/*!
\file strudyn.H
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

#ifndef STRUDYN_H
#define STRUDYN_H


#ifdef __cplusplus /* stuff which is c++ specific */

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

using namespace std;
using namespace Teuchos;



#endif /* end of #ifdef __cplusplus */





/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C" 
{
#endif

void strudyn();

#ifdef __cplusplus
}
#endif






#endif  /* #ifndef STRUDYN_H */
#endif  /* #ifdef CCADISCRET */
