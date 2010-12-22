/*!----------------------------------------------------------------------
\file drt_dofset_transparent.cpp

\brief A set of degrees of freedom special for contact

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_dofset_transparent.H"

DRT::TransparentDofSet::TransparentDofSet(RCP<DRT::Discretization> sourcedis) :
DRT::DofSet(),
sourcedis_(sourcedis)
{
  return;
}

int DRT::TransparentDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start)
{

  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis,dspos,start);

  TransferDegreesOfFreedom(*sourcedis_, dis, start);

  // tell all proxies (again!)
  NotifyAssigned();

  return count;
}

/// Assign dof numbers for new discretization using dof numbering from source discretization.
void DRT::TransparentDofSet::TransferDegreesOfFreedom(
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

    //build dofrowmap
    vector<int> dofrowvec;
    dofrowvec.reserve(dofrowmap_->NumMyElements());
    for (int inode = 0; inode != newdis.NumMyRowNodes(); ++inode)
    {
        const DRT::Node* newnode = newdis.lRowNode(inode);
        const DRT::Node* sourcenode = sourcedis.gNode(newnode->Id());

        const vector<int> dofs = sourcedis.Dof(sourcenode);

        const int newlid = newnode->LID();
        const int numdofs = (*numdfcolnodes_)[newlid];
        if (numdofs > 0)
        {
          (*idxcolnodes_)[newlid] = dofs[0];
          for (int idof = 0; idof < numdofs; ++idof)
          {
            dofrowvec.push_back(dofs[idof]);
          }
        }
    }
    dofrowmap_ = rcp(new Epetra_Map(-1, dofrowvec.size(), &dofrowvec[0], 0, newdis.Comm()));

    //build dofcolvec
    vector<int> dofcolvec;
    dofcolvec.reserve(dofcolmap_->NumMyElements());
    for (int inode = 0; inode != newdis.NumMyColNodes(); ++inode)
    {
        const DRT::Node* newnode = newdis.lColNode(inode);
        const DRT::Node* sourcenode = sourcedis.gNode(newnode->Id());
        const vector<int> dofs = sourcedis.Dof(sourcenode);
        const int newlid = newnode->LID();
        //const int newfirstidx = (*idxcolnodes_)[newlid];
        const int numdofs = (*numdfcolnodes_)[newlid];
        if (numdofs > 0)
        {
          (*idxcolnodes_)[newlid] = dofs[0];
          for (int idof = 0; idof < numdofs; ++idof)
          {
            dofcolvec.push_back(dofs[idof]);
          }
        }
    }
    dofcolmap_ = rcp(new Epetra_Map(-1, dofcolvec.size(), &dofcolvec[0], 0, newdis.Comm()));
}

#endif  // #ifdef CCADISCRET
