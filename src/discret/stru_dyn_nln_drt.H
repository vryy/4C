/*!----------------------------------------------------------------------
\file global_inp_control2.H
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

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
#include "drt_discret.H"
#include "linalg_utils.H"
#include "linalg_solver.H"




extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

#include "dstrc.H" /* needs standardtypes.h */

using namespace std;
using namespace Teuchos;



#endif /* end of #ifdef __cplusplus */





/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C" 
{
#endif

void dyn_nlnstructural_drt();
void stru_genalpha_drt();

#ifdef __cplusplus
}
#endif






#endif  /* #ifndef STRU_DYN_NLN_DRT_H */
#endif  /* #ifdef TRILINOS_PACKAGE */
#endif  /* #ifdef CCADISCRET */
