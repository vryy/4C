/*!----------------------------------------------------------------------
\file mortar_dofset.cpp
\brief A set of degrees of freedom special for mortar coupling

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

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
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

//#include <iostream>
//#include <algorithm>
//#include <numeric>

#include "mortar_dofset.H"
#include "mortar_node.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
MORTAR::MortarDofSet::MortarDofSet() :
DRT::DofSet()
{
  return;
}

/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int MORTAR::MortarDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis,
                                          const unsigned dspos, const int start)
{
  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis,dspos,start);

  // we'll get ourselves the row and column dof maps from the base class
  // and later replace them with our own version of them
  int nummyrow = dofrowmap_->NumMyElements();
  vector<int> myrow(nummyrow);
  int nummycol = dofcolmap_->NumMyElements();
  vector<int> mycol(nummycol);

  // now we loop all nodes in dis and create the new dof vectors
  for (int i=0; i<dis.NumMyColNodes(); ++i)
  {
    DRT::Node* node = dis.lColNode(i);
    if (!node) dserror("Cannot find local column node %d",i);
    // get dofs of node as created by base class DofSet
    vector<int> gdofs = Dof(node);
    // get dofs of node as we want them
    MORTAR::MortarNode* mrtrnode = dynamic_cast<MORTAR::MortarNode*>(node);
    if (!mrtrnode) dserror("dynamic_cast DRT::Node -> MORTAR::MortarNode failed");
    const int* newdofs = mrtrnode->Dofs();
    for (int j=0; j<(int)gdofs.size(); ++j)
    {
      if (!dofcolmap_->MyGID(gdofs[j])) dserror("Mismatch in degrees of freedom");
      int lid = dofcolmap_->LID(gdofs[j]);
      mycol[lid] = newdofs[j];
      // build dofrowmap as well
      if (!dofrowmap_->MyGID(gdofs[j])) continue;
      lid = dofrowmap_->LID(gdofs[j]);
      myrow[lid] = newdofs[j];
    }
  }

  // we have new vectors, so recreate epetra maps and replace old ones with them
  RCP<Epetra_Map> newdofrowmap = rcp(new Epetra_Map(-1,nummyrow,&myrow[0],0,dofrowmap_->Comm()));
  RCP<Epetra_Map> newdofcolmap = rcp(new Epetra_Map(-1,nummycol,&mycol[0],0,dofcolmap_->Comm()));

  // be a little psychotic in checking whether everything is ok....
  if (newdofrowmap->NumMyElements() != dofrowmap_->NumMyElements() ||
      newdofrowmap->NumGlobalElements() != dofrowmap_->NumGlobalElements() ||
      newdofcolmap->NumMyElements() != dofcolmap_->NumMyElements() ||
      newdofcolmap->NumGlobalElements() != dofcolmap_->NumGlobalElements() ||
      !newdofrowmap->UniqueGIDs())
    dserror("Something's wrong in dof maps");

  // replace the old maps by our new ones (note: automatically deletes old ones)
  dofrowmap_ = newdofrowmap;
  dofcolmap_ = newdofcolmap;

  // tell all proxies (again!)
  NotifyAssigned();

  return count;
}

#endif  // #ifdef CCADISCRET
