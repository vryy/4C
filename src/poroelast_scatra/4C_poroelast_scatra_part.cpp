// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_poroelast_scatra_part.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
PoroElastScaTra::PoroScatraPart::PoroScatraPart(
    MPI_Comm comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraBase(comm, timeparams)
{
  poro_field()->setup_solver();
}

FOUR_C_NAMESPACE_CLOSE
