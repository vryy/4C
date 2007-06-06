/*!----------------------------------------------------------------------
\file condif_drt.H
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifdef D_FLUID

#ifndef CONDIF_DRT_H
#define CONDIF_DRT_H



/* BEGIN   C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/
#ifdef __cplusplus /* stuff which is c++ specific */

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"
#include "../io/io_drt.H"

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
}

#include "../drt_lib/dstrc.H" /* needs standardtypes.h */

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
void dyn_condif_drt();



    
/* BEGIN   C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/
#ifdef __cplusplus
}
#endif
/* END     C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/


#endif  /* #ifndef CONDIF_DRT_H */
#endif  /* #ifdef D_FLUID */
#endif  /* #ifdef TRILINOS_PACKAGE */
#endif  /* #ifdef CCADISCRET */
