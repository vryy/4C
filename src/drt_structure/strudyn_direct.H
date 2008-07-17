/*----------------------------------------------------------------------*/
/*!
\file strudyn_direct.H
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

#ifndef STRUDYN_DIRECT_H
#define STRUDYN_DIRECT_H


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

#include "strutimint.H"
#include "strutimint_impl.H"
#include "strutimint_expl.H"
#include "strutimint_genalpha.H"
#include "strutimint_ost.H"
#include "strutimint_ab2.H"

#include "strutimada.H"

using namespace std;
using namespace Teuchos;



#endif /* end of #ifdef __cplusplus */





/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C" 
{
#endif

/*----------------------------------------------------------------------*/
//! Direct time integration of temporally (and spatially) discretised
//! structure
//!
//! \author bborn \date 06/08
void strudyn_direct();

/*----------------------------------------------------------------------*/
//! Create marching time integrator convenience routine
//!
//! \author bborn \date 07/08 
Teuchos::RCP<STR::StruTimInt> strudyn_CreateMarching
(
  const Teuchos::ParameterList& ioflags,  //!< input-output-flags
  const Teuchos::ParameterList& sdyn,  //!< structural dynamic flags
  const Teuchos::ParameterList& xparams,  //!< extra flags
  RefCountPtr<DRT::Discretization>& actdis,  //!< discretisation
  LINALG::Solver&  solver,  //!< the solver
  IO::DiscretizationWriter& output  //!< output writer
);

/*----------------------------------------------------------------------*/
//! Create auxiliar time integrator convenience routine
//!
//! \author bborn \date 07/08
Teuchos::RCP<STR::StruTimAda> strudyn_CreateAuxiliar
(
  const Teuchos::ParameterList& sdyn,  //!< TIS input parameters
  const Teuchos::ParameterList& tap,  //!< adaptive input flags
  Teuchos::RCP<STR::StruTimInt> tis  //!< marching time integrator
);

#ifdef __cplusplus
}
#endif



/*----------------------------------------------------------------------*/
#endif  /* #ifndef STRUDYN_DIRECT_H */
#endif  /* #ifdef CCADISCRET */
