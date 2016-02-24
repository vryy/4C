/*!----------------------------------------------------------------------
\file  drt_potential_dofset.cpp

\brief A set of degrees of freedom for potential discretizations

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#include "drt_potential_dofset.H"
#include "../linalg/linalg_utils.H"


/// constructor
POTENTIAL::PotentialDofSet::PotentialDofSet (Teuchos::RCP<DRT::Discretization> sourcedis) :
DRT::DofSet(),
sourcedis_(sourcedis)
{
  return;
}


/// Assign dof numbers for new discretization using dof numbering from source discretization.
int POTENTIAL::PotentialDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start)
{

  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis,dspos,start);
  if (pccdofhandling_) dserror("ERROR: Point coupling cinditions not yet implemented for PotentialDofSet");

  TransferDegreesOfFreedom(*sourcedis_, dis, start);

  // tell all proxies (again!)
  NotifyAssigned();

  return count;
}

/// Transfer degrees of freedom
void POTENTIAL::PotentialDofSet::TransferDegreesOfFreedom(
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

    // build local dofrowmap from source discretization with identical ids
    int countrowdof = 0;
    Teuchos::RCP< Epetra_IntVector > localrowdofs = Teuchos::rcp(new Epetra_IntVector(*newdis.DofRowMap()));
    Teuchos::RCP< Epetra_IntVector > localcoldofs = Teuchos::rcp(new Epetra_IntVector(*newdis.DofColMap()));

    Epetra_IntVector idxrownodes(*newdis.NodeRowMap());

    for (int inode = 0; inode != newdis.NumMyRowNodes(); ++inode)
    {
      const DRT::Node* newnode = newdis.lRowNode(inode);
      if(!sourcedis.HaveGlobalNode(newnode->Id()))
        dserror("source dis does not have node");

      const DRT::Node* sourcenode = sourcedis.gNode(newnode->Id());
      const std::vector<int> dofs = sourcedis.Dof(sourcenode);

      if( (newnode->Owner() != sourcenode->Owner()) ||  (newnode->Owner() != newdis.Comm().MyPID()) )
        dserror("node not on proc");

      // returns local col map id
      const int newlid = newnode->LID();
      const int numdofs = (*numdfcolnodes_)[newlid];
      dsassert(sourcedis.NumDof(sourcenode)==newdis.NumDof(newnode), "number of dofs does not match!");
      std::copy(dofs.begin(),dofs.end(),&(*localrowdofs)[countrowdof]);
      countrowdof += numdofs;

      if (dofs.size()>0)
        idxrownodes[inode] = dofs[0];
    }

    Epetra_Import nodeimporter( idxcolnodes_->Map(), idxrownodes.Map() );
    int err = idxcolnodes_->Import( idxrownodes, nodeimporter, Insert );
    if (err) dserror( "Import using importer returned err=%d", err );

    // import localrowdofs into localcoldofs
    // in this way the the original dofcolmap is preserved
    Epetra_Import importer(localcoldofs->Map(),localrowdofs->Map());

    err = localcoldofs->Import((*localrowdofs),importer,Insert);
    if (err) dserror("Import using importer returned err=%d",err);

    dofrowmap_ = Teuchos::rcp(new Epetra_Map(-1,(*localrowdofs).MyLength(),&(*localrowdofs)[0],0,newdis.Comm()));
    if (!dofrowmap_->UniqueGIDs()) dserror("Dof row map is not unique");

    dofcolmap_ = Teuchos::rcp(new Epetra_Map(-1,(*localcoldofs).MyLength(),&(*localcoldofs)[0],0,newdis.Comm()));

}

