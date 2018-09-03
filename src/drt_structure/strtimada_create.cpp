/*----------------------------------------------------------------------*/
/*!
\file strtimada_create.cpp
\brief Creation of auxiliary time integration scheme for time step size adaptivity

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "strtimada_create.H"
#include "strtimada_zienxie.H"
#include "strtimada_joint.H"
#include "strtimint_ab2.H"
#include "strtimint_expleuler.H"
#include "strtimint_centrdiff.H"

/*======================================================================*/
/* create auxiliary time integration scheme */
Teuchos::RCP<STR::TimAda> STR::TimAdaCreate(
    const Teuchos::ParameterList& ioflags, const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& sdyn, const Teuchos::ParameterList& xparams,
    const Teuchos::ParameterList& tap,  //!< adaptive input flags
    Teuchos::RCP<STR::TimInt> tis       //!< marching time integrator
)
{
  Teuchos::RCP<STR::TimAda> sta = Teuchos::null;

  // auxiliary time integrator
  switch (DRT::INPUT::IntegralValue<INPAR::STR::TimAdaKind>(tap, "KIND"))
  {
    case INPAR::STR::timada_kind_none:
      // No adaptivity in time
      sta = Teuchos::null;
      break;

    case INPAR::STR::timada_kind_zienxie:
      // Zienkiewicz-Xie error indicator for generalised-alpha
      sta = Teuchos::rcp(new STR::TimAdaZienXie(timeparams, tap, tis));
      break;

    case INPAR::STR::timada_kind_ab2:
      // Adams-Bashforth 2nd order
      sta = Teuchos::rcp(
          new STR::TimAdaJoint<STR::TimIntAB2>(ioflags, timeparams, sdyn, xparams, tap, tis));
      break;

    case INPAR::STR::timada_kind_expleuler:
      // Adams-Bashforth 2nd order
      sta = Teuchos::rcp(
          new STR::TimAdaJoint<STR::TimIntExplEuler>(ioflags, timeparams, sdyn, xparams, tap, tis));
      break;

    case INPAR::STR::timada_kind_centraldiff:
      // Adams-Bashforth 2nd order
      sta = Teuchos::rcp(
          new STR::TimAdaJoint<STR::TimIntCentrDiff>(ioflags, timeparams, sdyn, xparams, tap, tis));
      break;

    default:
      dserror("Auxiliary time integrator is not available.");
      break;
  }

  // return the auxiliary integrator
  return sta;
}

/*----------------------------------------------------------------------*/
