/*!----------------------------------------------------------------------
\file drt_cdofset.cpp

\brief A set of degrees of freedom special for contact

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

//#include <iostream>
//#include <algorithm>
//#include <numeric>

#include "drt_cdofset.H"
#include "drt_cnode.H"
/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
CONTACT::CDofSet::CDofSet() :
DRT::DofSet()
{
  return;
}


/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int CONTACT::CDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const int start)
{
  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis,start);
  
  // now we'll get ourselves the row and column dof maps from the base class and replace
  // them with our own version of them
  //RCP<Epetra_Map> oldrowmap = dofrowmap_;
  //RCP<Epetra_Map> oldcolmap = dofcolmap_;

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
    CONTACT::CNode* cnode = dynamic_cast<CONTACT::CNode*>(node);
    if (!cnode) dserror("dynamic_cast DRT::Node -> CONTACT::CNode failed");
    const int* newdofs = cnode->Dofs();
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
    
  return count;
}


#endif  // #ifdef CCADISCRET
