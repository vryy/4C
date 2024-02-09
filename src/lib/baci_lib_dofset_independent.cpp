/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of independent dofset

\level 2


*/
/*---------------------------------------------------------------------*/

#include "baci_lib_dofset_independent.hpp"

#include "baci_lib_discret.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::IndependentDofSet::IndependentDofSet(bool ignoreminnodegid /*=false*/)
    : DRT::DofSet(), ignoreminnodegid_(ignoreminnodegid)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::IndependentDofSet::IndependentDofSet(const IndependentDofSet& old)
    : DRT::DofSet(old), ignoreminnodegid_(old.ignoreminnodegid_)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::IndependentDofSet::AddDofSettoList()
{
  // We do nothing here as an independent DofSet should not show up in the dof set list.
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::IndependentDofSet::GetFirstGIDNumberToBeUsed(const Discretization& dis) const
{
  // always start from zero here, as this is an independent DofSet
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::IndependentDofSet::GetMinimalNodeGIDIfRelevant(const Discretization& dis) const
{
  return ignoreminnodegid_ ? 0 : dis.NodeRowMap()->MinAllGID();
}

BACI_NAMESPACE_CLOSE
