// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_dofset_independent.hpp"

#include "4C_fem_discretization.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::IndependentDofSet::IndependentDofSet(bool ignoreminnodegid /*=false*/)
    : Core::DOFSets::DofSet(), ignoreminnodegid_(ignoreminnodegid)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::IndependentDofSet::IndependentDofSet(const IndependentDofSet& old)
    : Core::DOFSets::DofSet(old), ignoreminnodegid_(old.ignoreminnodegid_)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::IndependentDofSet::add_dof_setto_list()
{
  // We do nothing here as an independent DofSet should not show up in the dof set list.
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::DOFSets::IndependentDofSet::get_first_gid_number_to_be_used(
    const Core::FE::Discretization& dis) const
{
  // always start from zero here, as this is an independent DofSet
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::DOFSets::IndependentDofSet::get_minimal_node_gid_if_relevant(
    const Core::FE::Discretization& dis) const
{
  return ignoreminnodegid_ ? 0 : dis.node_row_map()->MinAllGID();
}

FOUR_C_NAMESPACE_CLOSE
