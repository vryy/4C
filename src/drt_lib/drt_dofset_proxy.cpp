/*----------------------------------------------------------------------*/
/*!
\file drt_dofset_proxy.cpp

\brief Implementation Proxy to a set of degrees of freedom

<pre>
\brief Implementation
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#include "drt_dofset_proxy.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetProxy::DofSetProxy(DofSet* dofset)
  : dofset_(dofset)
{
  dofset->RegisterProxy(this);
  NotifyAssigned();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetProxy::~DofSetProxy()
{
  if (dofset_!=NULL)
    dofset_->UnregisterProxy(this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetProxy::AddDofSettoList()
{
  // We do nothing here as a proxy does not show up in the dof set list.
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetProxy::AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start)
{
  // Assume our original DofSet is valid right now. Otherwise we will be
  // notified anyway.
  NotifyAssigned();
  return start;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetProxy::NotifyAssigned()
{
  if (dofset_!=NULL)
  {
    // Just copy those Teuchos::rcps.
    dofrowmap_        = dofset_->dofrowmap_;
    dofcolmap_        = dofset_->dofcolmap_;
    numdfcolnodes_    = dofset_->numdfcolnodes_;
    numdfcolelements_ = dofset_->numdfcolelements_;
    idxcolnodes_      = dofset_->idxcolnodes_;
    idxcolelements_   = dofset_->idxcolelements_;

    // special dof handling due to point coupling conditions
    pccdofhandling_   = dofset_->pccdofhandling_;
    dofscolnodes_     = dofset_->dofscolnodes_;
    shiftcolnodes_    = dofset_->shiftcolnodes_;
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetProxy::NotifyReset()
{
  // clear my Teuchos::rcps.
  Reset();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetProxy::Disconnect(DofSet* dofset)
{
  if (dofset==dofset_)
    dofset_ = NULL;
  else
    dserror("cannot disconnect from non-connected DofSet");

  // clear my Teuchos::rcps.
  Reset();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::DofSetProxy::Filled() const
{
  if (dofset_)
  {
    return dofset_->Filled();
  }
  return false;
}

