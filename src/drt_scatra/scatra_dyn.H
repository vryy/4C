/*!----------------------------------------------------------------------
\file scatra_dyn.H
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef SCATRA_DYN_H
#define SCATRA_DYN_H


/* BEGIN   C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/
#ifdef __cplusplus /* stuff which is c++ specific */

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */


using namespace std;
using namespace Teuchos;


#endif /* end of #ifdef __cplusplus */
/* END     C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/


/* BEGIN   C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/
#ifdef __cplusplus
extern "C"
{
#endif
/* END     C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/


/* headers and prototypes of functions also acessible from c */

/*! entry point for the solution of scalar transport problems */
void scatra_dyn(
    int disnumff,      /* fluid discretization number, typically 0 */
    int disnumscatra,  /* scatra discretization number, typically 1 */
    int restart        /* do we have to perform a restart?  */
    );



/* BEGIN   C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/
#ifdef __cplusplus
}
#endif
/* END     C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/


#endif  /* #ifndef SCATRA_DYN_H */
#endif  /* #ifdef CCADISCRET */
