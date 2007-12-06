/*!----------------------------------------------------------------------
\file xfluid_dyn_nln_drt.H
\brief Control routine for fluid (time) integration. Calls

     o Single step one-step-theta time integration

     o Two step BDF2 Gear's methode with one-step-theta start step

     o Genalpha time integartion

     o Stationary solver


<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef XFLUID_DYN_NLN_DRT_H
#define XFLUID_DYN_NLN_DRT_H



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
void xdyn_fluid_drt();



    
/* BEGIN   C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/
#ifdef __cplusplus
}
#endif
/* END     C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/


#endif  /* #ifndef XFLUID_DYN_NLN_DRT_H */
#endif  /* #ifdef CCADISCRET */
