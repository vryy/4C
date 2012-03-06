/*!----------------------------------------------------------------------

\file stru_dyn_nln_drt.H

\brief Control routine for structural dynamics (outsourced to adapter layer)

<pre>
Maintainer: Thomas Klöppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef STRU_DYN_NLN_DRT_H
#define STRU_DYN_NLN_DRT_H


#ifdef __cplusplus /* stuff which is c++ specific */

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"


using namespace std;
using namespace Teuchos;



#endif /* end of #ifdef __cplusplus */





/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C"
{
#endif

void dyn_nlnstructural_drt();
void caldyn_drt();

#ifdef __cplusplus
}
#endif






#endif  /* #ifndef STRU_DYN_NLN_DRT_H */
#endif  /* #ifdef CCADISCRET */
