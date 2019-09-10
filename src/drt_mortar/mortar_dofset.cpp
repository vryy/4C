/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of degrees of freedom special for mortar coupling

\maintainer Matthias Mayr

\level 1
*/
/*---------------------------------------------------------------------*/

#include "mortar_dofset.H"
#include "mortar_node.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
MORTAR::MortarDofSet::MortarDofSet() : DRT::DofSet() { return; }

/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int MORTAR::MortarDofSet::AssignDegreesOfFreedom(
    const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis, dspos, start);
  if (pccdofhandling_)
    dserror("ERROR: Point coupling cinditions not yet implemented for MortarDofSet");

  // we'll get ourselves the row and column dof maps from the base class
  // and later replace them with our own version of them
  int nummyrow = dofrowmap_->NumMyElements();
  std::vector<int> myrow(nummyrow);
  int nummycol = dofcolmap_->NumMyElements();
  std::vector<int> mycol(nummycol);

  // now we loop all nodes in dis and create the new dof vectors
  for (int i = 0; i < dis.NumMyColNodes(); ++i)
  {
    DRT::Node* node = dis.lColNode(i);
    if (!node) dserror("Cannot find local column node %d", i);
    // get dofs of node as created by base class DofSet
    std::vector<int> gdofs = Dof(node);
    // get dofs of node as we want them
    MORTAR::MortarNode* mrtrnode = dynamic_cast<MORTAR::MortarNode*>(node);
    if (!mrtrnode) dserror("dynamic_cast DRT::Node -> MORTAR::MortarNode failed");
    const int* newdofs = mrtrnode->Dofs();
    for (int j = 0; j < (int)gdofs.size(); ++j)
    {
      if (!dofcolmap_->MyGID(gdofs[j])) dserror("Mismatch in degrees of freedom");
      int lid = dofcolmap_->LID(gdofs[j]);
      mycol[lid] = newdofs[j];
      // build dofrowmap as well
      if (!dofrowmap_->MyGID(gdofs[j])) continue;
      lid = dofrowmap_->LID(gdofs[j]);
      myrow[lid] = newdofs[j];
    }
    if (gdofs.size() > 0)
    {
      (*idxcolnodes_)[i] = newdofs[0];
    }
  }

  // we have new vectors, so recreate epetra maps and replace old ones with them
  Teuchos::RCP<Epetra_Map> newdofrowmap =
      Teuchos::rcp(new Epetra_Map(-1, nummyrow, &myrow[0], 0, dofrowmap_->Comm()));
  Teuchos::RCP<Epetra_Map> newdofcolmap =
      Teuchos::rcp(new Epetra_Map(-1, nummycol, &mycol[0], 0, dofcolmap_->Comm()));

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

  // mortar element == face element: we need this...
  idxcolfaces_ = idxcolelements_;
  numdfcolfaces_ = numdfcolelements_;

  // tell all proxies (again!)
  NotifyAssigned();

  return count;
}
