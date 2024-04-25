/*----------------------------------------------------------------------*/
/*! \file
\brief Creation of auxiliary time integration scheme for time step size adaptivity
\level 1
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_structure_timada_create.hpp"

#include "4C_structure_timada_joint.hpp"
#include "4C_structure_timada_zienxie.hpp"
#include "4C_structure_timint_ab2.hpp"
#include "4C_structure_timint_centrdiff.hpp"
#include "4C_structure_timint_expleuler.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>

FOUR_C_NAMESPACE_OPEN

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
  switch (CORE::UTILS::IntegralValue<INPAR::STR::TimAdaKind>(tap, "KIND"))
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
      FOUR_C_THROW("Auxiliary time integrator is not available.");
      break;
  }

  // return the auxiliary integrator
  return sta;
}

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
