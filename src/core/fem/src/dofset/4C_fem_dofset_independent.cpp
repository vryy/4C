/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of independent dofset

\level 2


*/
/*---------------------------------------------------------------------*/

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
