/*----------------------------------------------------------------------*/
/*!
\file drt_dofset_proxy.cpp

\brief Proxy to a set of degrees of freedom

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
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
  }
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

