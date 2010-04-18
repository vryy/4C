/*!----------------------------------------------------------------------
\file drt_dofset_base.cpp
\brief A set of degrees of freedom

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
Maintainer: Ulrrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <iostream>
#include <algorithm>
#include <numeric>

#include "drt_dofset_base.H"
#include "drt_discret.H"
#include "drt_utils.H"

#include "linalg_utils.H"

// list of all dof sets
std::list<DRT::DofSetBase*> DRT::DofSetBase::static_dofsets_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSetBase::DofSetBase()
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSetBase::~DofSetBase()
{
  static_dofsets_.remove(this);
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetBase::NumGlobalElements() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSetBase::NumGlobalElements(): dofrowmap_ not initialized, yet");
  return dofrowmap_->NumGlobalElements();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetBase::MaxAllGID() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSetBase::MaxAllGID(): dofrowmap_ not initialized, yet");
  return dofrowmap_->MaxAllGID();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DofSetBase::DofRowMap() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSetBase::DofRowMap(): dofrowmap_ not initialized, yet");
  return dofrowmap_.get();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DofSetBase::DofColMap() const
{
  if (dofcolmap_ == Teuchos::null)
    dserror("DRT::DofSetBase::DofColMap(): dofcolmap_ not initialized, yet");
  return dofcolmap_.get();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetBase::AddDofSettoList()
{
  if (std::find(static_dofsets_.begin(),static_dofsets_.end(),this)==static_dofsets_.end())
  {
    static_dofsets_.push_back(this);
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::DofSetBase::Initialized() const
{
  if (dofcolmap_ == Teuchos::null or dofrowmap_ == Teuchos::null)
    return false;
  else
    return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetBase::MaxGIDinList() const
{
  int count = -1;
  for (std::list<DofSetBase*>::const_iterator i=static_dofsets_.begin();
       i!=static_dofsets_.end();
       ++i)
  {
    if (*i==this)
      break;

    // ignore empty (no yet initialized) dof row maps
    if ((*i)->NumGlobalElements()>0)
    {
      count = max((*i)->MaxAllGID(),count);
    }
  }
  return count;
}

#endif  // #ifdef CCADISCRET
