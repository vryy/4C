/*======================================================================*/
/*!
\file stru_genalpha_zienxie_drt.H

\brief Structural time integration with generalised-alpha (marching scheme)
       and Zienkiewicz-Xie (auxiliar scheme) as error indicator.

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef STRU_GENALPHA_ZIENXIE_H
#define STRU_GENALPHA_ZIENXIE_H


/*----------------------------------------------------------------------*/
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




extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

using namespace std;
using namespace Teuchos;

   void genalpha_setparalist
   (
      const bool damp,
      const double damp_k,
      const double damp_m,
      const double beta,
      const double gamma,
      const double alpham,
      const double alphaf,
      const double timein,
      const double stepsize,
      const int step,
      const int stepmax,
      const int itermax,
      const int iter,
      const double toldisp,
      const bool ioflags_struct_disp,
      const int updevry_disp,
      const bool ioflags_struct_stress,
      const int updevry_stress,
      const int restart,
      const int res_write_evry,
      const bool print_to_stdout,
      const bool print_to_errfile,
      FILE* errfilein,
      const STRUCT_DYNAMIC::_nlnSolvTyp nst,
      const string pred,
      ParameterList& params
   );

#endif /* end of #ifdef __cplusplus */



/*----------------------------------------------------------------------*/
/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C" 
{
#endif

   void stru_genalpha_zienxie_drt();

#ifdef __cplusplus
}
#endif






#endif  /* #ifndef STRU_GENALPHA_ZIENXIE_H */
#endif  /* #ifdef CCADISCRET */
