/*----------------------------------------------------------------------*/
/*!
\file strtimada_create.H
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

#ifndef STRTIMADA_CREATE_H
#define STRTIMADA_CREATE_H

/*----------------------------------------------------------------------*/
/* headers */
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

#include "strtimint.H"
#include "strtimint_impl.H"
#include "strtimint_expl.H"
#include "strtimint_genalpha.H"
#include "strtimint_ost.H"
#include "strtimint_gemm.H"
#include "strtimint_ab2.H"

#include "strtimada.H"
#include "strtimada_zienxie.H"

/*----------------------------------------------------------------------*/
/* Structural Dynamics */
namespace STR
{

  /*====================================================================*/
  //! Create auxiliar time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimAda> TimAdaCreate
  (
    const Teuchos::ParameterList& ioflags,  //!< input-output-flags
    const Teuchos::ParameterList& sdyn,  //!< structural dynamic flags
    const Teuchos::ParameterList& xparams,  //!< extra flags
    const Teuchos::ParameterList& tap,  //!< adaptive input flags
    Teuchos::RCP<STR::TimInt> tis  //!< marching time integrator
  );

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  /* #ifndef STRTIMADA_CREATE_H */
#endif  /* #ifdef CCADISCRET */
