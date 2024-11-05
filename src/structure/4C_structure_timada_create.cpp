// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_timada_create.hpp"

#include "4C_structure_timada_joint.hpp"
#include "4C_structure_timada_zienxie.hpp"
#include "4C_structure_timint_ab2.hpp"
#include "4C_structure_timint_centrdiff.hpp"
#include "4C_structure_timint_expleuler.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* create auxiliary time integration scheme */
std::shared_ptr<Solid::TimAda> Solid::tim_ada_create(
    const Teuchos::ParameterList& ioflags, const Teuchos::ParameterList& timeparams,
    const Teuchos::ParameterList& sdyn, const Teuchos::ParameterList& xparams,
    const Teuchos::ParameterList& tap,  //!< adaptive input flags
    std::shared_ptr<Solid::TimInt> tis  //!< marching time integrator
)
{
  std::shared_ptr<Solid::TimAda> sta = nullptr;

  // auxiliary time integrator
  switch (Teuchos::getIntegralValue<Inpar::Solid::TimAdaKind>(tap, "KIND"))
  {
    case Inpar::Solid::timada_kind_none:
      // No adaptivity in time
      sta = nullptr;
      break;

    case Inpar::Solid::timada_kind_zienxie:
      // Zienkiewicz-Xie error indicator for generalised-alpha
      sta = std::make_shared<Solid::TimAdaZienXie>(timeparams, tap, tis);
      break;

    case Inpar::Solid::timada_kind_ab2:
      // Adams-Bashforth 2nd order
      sta = std::make_shared<Solid::TimAdaJoint<Solid::TimIntAB2>>(
          ioflags, timeparams, sdyn, xparams, tap, tis);
      break;

    case Inpar::Solid::timada_kind_expleuler:
      // Adams-Bashforth 2nd order
      sta = std::make_shared<Solid::TimAdaJoint<Solid::TimIntExplEuler>>(
          ioflags, timeparams, sdyn, xparams, tap, tis);
      break;

    case Inpar::Solid::timada_kind_centraldiff:
      // Adams-Bashforth 2nd order
      sta = std::make_shared<Solid::TimAdaJoint<Solid::TimIntCentrDiff>>(
          ioflags, timeparams, sdyn, xparams, tap, tis);
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
