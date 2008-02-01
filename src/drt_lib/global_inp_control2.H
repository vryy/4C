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

#ifndef GLOBAL_INP_CONTROL2_H
#define GLOBAL_INP_CONTROL2_H


#ifdef __cplusplus /* stuff which is c++ specific */

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"
#include "drt_inputreader.H"

#include "drt_globalproblem.H"
#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_condition.H"
#include "drt_utils.H"
#include "linalg_utils.H"
#ifdef D_SHELL8
#include "../drt_s8/shell8.H"
#endif
#ifdef D_FLUID2
#include "../drt_f2/fluid2.H"
#include "../drt_f2/condif2.H"
#endif
#ifdef D_FLUID3
#include "../drt_f3/fluid3.H"
#include "../drt_f3/xfluid3.H"
#endif


using namespace std;
using namespace Teuchos;

// headers and prototypes of functions acessible from c++

string inp_micro_filename(int mypid);

#endif /* end of #ifdef __cplusplus */



/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C"
{
#endif

void ntainp_ccadiscret();
void distribute_drt_grids();
void dyn_nlnstructural_drt();
void drt_problem_done();

#ifdef __cplusplus
}
#endif






#endif  /* #ifndef GLOBAL_INP_CONTROL2_H */

#endif  /* #ifdef CCADISCRET */
