/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of independent dofset

\level 2

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#include "drt_dofset_independent.H"

#include "../drt_lib/drt_discret.H"

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
DRT::IndependentDofSet::~IndependentDofSet() { return; }


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
