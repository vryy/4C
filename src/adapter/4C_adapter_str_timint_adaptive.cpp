// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_timint_adaptive.hpp"

#include "4C_structure_timada.hpp"
#include "4C_structure_timint.hpp"
#include "4C_structure_timint_create.hpp"

FOUR_C_NAMESPACE_OPEN


/*======================================================================*/
/* constructor */
Adapter::StructureTimIntAda::StructureTimIntAda(
    Teuchos::RCP<Solid::TimAda> sta, Teuchos::RCP<Structure> sti)
    : StructureWrapper(sti), sta_(sta)
{
  // make sure
  if (sta_ == Teuchos::null) FOUR_C_THROW("Failed to create structural integrator");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Adapter::StructureTimIntAda::integrate() { return sta_->integrate(); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimIntAda::prepare_output(bool force_prepare)
{
  sta_->prepare_output_period(force_prepare);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimIntAda::output() { sta_->output_period(); }


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
