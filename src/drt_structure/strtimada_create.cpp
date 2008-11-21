/*----------------------------------------------------------------------*/
/*!
\file strtimada_create.cpp
\brief Structural time integration

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

/*----------------------------------------------------------------------*/
/* headers */
#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "strudyn_direct.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "stru_resulttest.H"
#include "../drt_inv_analysis/inv_analysis.H"

#include "strtimint.H"
#include "strtimint_impl.H"
#include "strtimint_expl.H"
#include "strtimint_genalpha.H"
#include "strtimint_ost.H"
#include "strtimint_gemm.H"
#include "strtimint_ab2.H"

#include "strtimada.H"
#include "strtimada_zienxie.H"
#include "strtimada_joint.H"

#include "strtimada_create.H"


/*======================================================================*/
/* create auxiliar time integration scheme */
Teuchos::RCP<STR::TimAda> STR::TimAdaCreate
(
  const Teuchos::ParameterList& ioflags,
  const Teuchos::ParameterList& sdyn,
  const Teuchos::ParameterList& xparams,
  const Teuchos::ParameterList& tap,  //!< adaptive input flags
  Teuchos::RCP<STR::TimInt> tis  //!< marching time integrator
)
{
  Teuchos::RCP<STR::TimAda> sta = Teuchos::null;

  // auxiliar time integrator
  switch (Teuchos::getIntegralValue<int>(tap,"KIND"))
  {

  case TIMADA_DYNAMIC::timada_kind_none :
    // No adaptivity in time
    sta = Teuchos::null;
    break;

  case TIMADA_DYNAMIC::timada_kind_zienxie :
    // Zienkiewivz-Xie error indicator for generalised-alpha
    sta = Teuchos::rcp(new STR::TimAdaZienXie(sdyn, tap, tis));
    break;

  case TIMADA_DYNAMIC::timada_kind_ab2 :
    // Adams-Bashforth 2nd order
    sta = Teuchos::rcp(new STR::TimAdaJoint<STR::TimIntAB2>(ioflags, sdyn, xparams, tap, tis));
    break;

  default :
    dserror("Auxiliar time integrator is not available.");
    break;

  }

  // return the auxiliar integrator
  return sta;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
