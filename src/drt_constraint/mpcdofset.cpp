/*!----------------------------------------------------------------------
\file mpcdofset.cpp

\brief A set of degrees of freedom special for contact

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

//#include <iostream>
//#include <algorithm>
//#include <numeric>

#include "mpcdofset.H"
/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
MPCDofSet::MPCDofSet(RCP<DRT::Discretization> sourcedis) :
DRT::DofSet(),
sourcedis_(sourcedis)
{
  return;
}


/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int MPCDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const int start)
{
  // Important remark (popp, 04/08):
  // We explicitly set the flag "unique_and_unchanging_dofnumbers" = false!
  // This way we have no problems dealing with an arbitrary number of
  // contact interfaces. The idea for the flag = true comes from monolith.FSI
  // where we have multiple discretizations and want to have a unique
  // dof numbering across them. So, when we set the flag to false for contact
  // at the moment, this is no problem AS LONG AS we do not want to deal
  // with multiple discretizations! (FIXME)
  
  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis,start);
  
  TransferDegreesOfFreedom(*sourcedis_, dis, start);
    
  return count;
}

/// Assign dof numbers for new discretization using dof numbering from source discretization.
void MPCDofSet::TransferDegreesOfFreedom(
        const DRT::Discretization& sourcedis,
        const DRT::Discretization& newdis,
        const int start
        )
{
    if (!sourcedis.DofRowMap()->UniqueGIDs()) dserror("DofRowMap is not unique");
    if (!sourcedis.NodeRowMap()->UniqueGIDs()) dserror("NodeRowMap is not unique");
    if (!sourcedis.ElementRowMap()->UniqueGIDs()) dserror("ElementRowMap is not unique");

    if (!newdis.DofRowMap()->UniqueGIDs()) dserror("DofRowMap is not unique");
    if (!newdis.NodeRowMap()->UniqueGIDs()) dserror("NodeRowMap is not unique");
    if (!newdis.ElementRowMap()->UniqueGIDs()) dserror("ElementRowMap is not unique");
    
    
    if (!newdis.DofRowMap()->UniqueGIDs()) dserror("DofRowMap is not unique");
    if (!newdis.NodeRowMap()->UniqueGIDs()) dserror("NodeRowMap is not unique");
    if (!newdis.ElementRowMap()->UniqueGIDs()) dserror("ElementRowMap is not unique");
    
    if (!sourcedis.DofRowMap()->UniqueGIDs()) dserror("DofRowMap is not unique");
    if (!sourcedis.NodeRowMap()->UniqueGIDs()) dserror("NodeRowMap is not unique");
    if (!sourcedis.ElementRowMap()->UniqueGIDs()) dserror("ElementRowMap is not unique");

    if (!newdis.DofRowMap()->UniqueGIDs()) dserror("DofRowMap is not unique");
    if (!newdis.NodeRowMap()->UniqueGIDs()) dserror("NodeRowMap is not unique");
    if (!newdis.ElementRowMap()->UniqueGIDs()) dserror("ElementRowMap is not unique");

    
    vector<int> dofrowvec(dofrowmap_->NumMyElements(),1000);
    // now for the nodes
    int counter = 0;
    for (int inode = 0; inode != newdis.NumMyRowNodes(); ++inode)
    {
        const DRT::Node* newnode = newdis.lRowNode(inode);
        const DRT::Node* sourcenode = sourcedis.gNode(newnode->Id());
        
        const vector<int> dofs = sourcedis.Dof(sourcenode);
        dsassert(dofs.size()==3, "number of dofs is not 3!");
        
        const int newlid = newnode->LID();
        const int numdofs = (*numdfcolnodes_)[newlid];
        dsassert(numdofs==3, "number of dofs is not 3!");
        for (int idof = 0; idof < numdofs; ++idof)
        {
            dofrowvec[counter] = dofs[idof];
            counter++;
        }
    }
  
    
    if (!dofrowmap_->UniqueGIDs()) dserror("before, Dof row map is not unique");
    dofrowmap_ = rcp(new Epetra_Map(-1, dofrowvec.size(), &dofrowvec[0], 0, newdis.Comm()));
    if (!dofrowmap_->UniqueGIDs()) dserror("after, Dof row map is not unique");
    
    vector<int> dofcolvec(dofcolmap_->NumMyElements());
    int colcounter = 0;
    for (int inode = 0; inode != newdis.NumMyColNodes(); ++inode)
    {
        const DRT::Node* newnode = newdis.lColNode(inode);
        const DRT::Node* sourcenode = sourcedis.gNode(newnode->Id());
        const vector<int> dofs = sourcedis.Dof(sourcenode);
        const int newlid = newnode->LID();
        //const int newfirstidx = (*idxcolnodes_)[newlid];
        const int numdofs = (*numdfcolnodes_)[newlid];
        for (int idof = 0; idof < numdofs; ++idof)
        {
            dofcolvec[colcounter] = dofs[idof];
            colcounter++;
        }
    }
    dofcolmap_ = rcp(new Epetra_Map(-1, dofcolvec.size(), &dofcolvec[0], 0, newdis.Comm()));
}

#endif  // #ifdef CCADISCRET
