/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of independent dofset

\level 2


*/
/*---------------------------------------------------------------------*/

#include "4C_discretization_dofset_independent.hpp"

#include "4C_lib_discret.hpp"

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
void Core::DOFSets::IndependentDofSet::AddDofSettoList()
{
  // We do nothing here as an independent DofSet should not show up in the dof set list.
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::DOFSets::IndependentDofSet::get_first_gid_number_to_be_used(
    const Discret::Discretization& dis) const
{
  // always start from zero here, as this is an independent DofSet
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::DOFSets::IndependentDofSet::get_minimal_node_gid_if_relevant(
    const Discret::Discretization& dis) const
{
  return ignoreminnodegid_ ? 0 : dis.NodeRowMap()->MinAllGID();
}

FOUR_C_NAMESPACE_CLOSE
