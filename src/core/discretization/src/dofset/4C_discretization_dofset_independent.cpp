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
CORE::Dofsets::IndependentDofSet::IndependentDofSet(bool ignoreminnodegid /*=false*/)
    : CORE::Dofsets::DofSet(), ignoreminnodegid_(ignoreminnodegid)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Dofsets::IndependentDofSet::IndependentDofSet(const IndependentDofSet& old)
    : CORE::Dofsets::DofSet(old), ignoreminnodegid_(old.ignoreminnodegid_)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Dofsets::IndependentDofSet::AddDofSettoList()
{
  // We do nothing here as an independent DofSet should not show up in the dof set list.
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::Dofsets::IndependentDofSet::get_first_gid_number_to_be_used(
    const DRT::Discretization& dis) const
{
  // always start from zero here, as this is an independent DofSet
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::Dofsets::IndependentDofSet::get_minimal_node_gid_if_relevant(
    const DRT::Discretization& dis) const
{
  return ignoreminnodegid_ ? 0 : dis.NodeRowMap()->MinAllGID();
}

FOUR_C_NAMESPACE_CLOSE
