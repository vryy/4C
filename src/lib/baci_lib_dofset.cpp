/*---------------------------------------------------------------------*/
/*! \file

\brief A set of degrees of freedom

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_lib_dofset.H"

#include "baci_comm_exporter.H"
#include "baci_lib_discret.H"
#include "baci_lib_discret_hdg.H"
#include "baci_linalg_utils_sparse_algebra_math.H"

#include <Epetra_FECrsGraph.h>

#include <algorithm>
#include <iostream>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSet::DofSet() : DRT::DofSetBase(), filled_(false), dspos_(0), pccdofhandling_(false)
{
  return;
}



/*----------------------------------------------------------------------*
 |  << operator                                               ukue 04/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const DRT::DofSet& dofset)
{
  dofset.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this  (public)                                      ukue 04/07|
 *----------------------------------------------------------------------*/
void DRT::DofSet::Print(std::ostream& os) const
{
  for (int proc = 0; proc < numdfcolelements_->Comm().NumProc(); ++proc)
  {
    if (proc == numdfcolelements_->Comm().MyPID())
    {
      if (numdfcolelements_->MyLength()) os << "-------------------------- Proc " << proc << " :\n";
      for (int i = 0; i < numdfcolelements_->MyLength(); ++i)
      {
        int numdf = (*numdfcolelements_)[i];
        int idx = (*idxcolelements_)[i];
        os << i << ": ";
        for (int j = 0; j < numdf; ++j) os << (idx + j) << " ";
        os << "\n";
      }
      os << std::endl;
    }
    numdfcolelements_->Comm().Barrier();
  }
  for (int proc = 0; proc < numdfcolnodes_->Comm().NumProc(); ++proc)
  {
    if (proc == numdfcolnodes_->Comm().MyPID())
    {
      if (numdfcolnodes_->MyLength()) os << "-------------------------- Proc " << proc << " :\n";
      for (int i = 0; i < numdfcolnodes_->MyLength(); ++i)
      {
        int numdf = (*numdfcolnodes_)[i];
        int idx = (*idxcolnodes_)[i];
        os << i << ": ";
        for (int j = 0; j < numdf; ++j) os << (idx + j) << " ";
        os << "\n";
      }
      os << std::endl;
    }
    numdfcolnodes_->Comm().Barrier();
  }
  for (int proc = 0; proc < numdfcolfaces_->Comm().NumProc(); ++proc)
  {
    if (proc == numdfcolfaces_->Comm().MyPID())
    {
      if (numdfcolfaces_->MyLength()) os << "-------------------------- Proc " << proc << " :\n";
      for (int i = 0; i < numdfcolfaces_->MyLength(); ++i)
      {
        int numdf = (*numdfcolfaces_)[i];
        int idx = (*idxcolfaces_)[i];
        os << i << ": ";
        for (int j = 0; j < numdf; ++j) os << (idx + j) << " ";
        os << "\n";
      }
      os << std::endl;
    }
    numdfcolfaces_->Comm().Barrier();
  }
}


/*----------------------------------------------------------------------*
 |  reset everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
void DRT::DofSet::Reset()
{
  dofrowmap_ = Teuchos::null;
  dofcolmap_ = Teuchos::null;
  numdfcolnodes_ = Teuchos::null;
  numdfcolelements_ = Teuchos::null;
  idxcolnodes_ = Teuchos::null;
  idxcolelements_ = Teuchos::null;
  shiftcolnodes_ = Teuchos::null;
  dofscolnodes_ = Teuchos::null;

  filled_ = false;

  // tell all proxies
  NotifyReset();
}

/*----------------------------------------------------------------------*
 |  setup everything  (public)                                ukue 04/07|
 *----------------------------------------------------------------------*/
int DRT::DofSet::AssignDegreesOfFreedom(
    const Discretization& dis, const unsigned dspos, const int start)
{
  if (!dis.Filled()) dserror("discretization Filled()==false");
  if (!dis.NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!dis.ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

  // A definite offset is currently not supported.
  // TODO (kronbichler) find a better solution for this
  // if (start!=0)
  //  dserror("right now user specified dof offsets are not supported");

  dspos_ = dspos;

  // Add DofSets in order of assignment to list. Once it is there it has its
  // place and will get its starting id from the previous DofSet.
  AddDofSettoList();

  // We assume that all dof sets before this one have been set up. Otherwise
  // we'd have to reorder the list.
  //
  // There is no test anymore to make sure that all prior dof sets have been
  // assigned. It seems people like to manipulate dof sets. People do create
  // dof sets that do not contain any dofs (on its first assignment), people
  // even shift dof set numbers to create overlapping dof sets. This is
  // perfectly fine.
  //
  // However if you rely on non-overlapping dof sets, you have to
  // FillComplete() your discretizations in the order of their creation. This
  // is guaranteed for all discretizations read from the input file since the
  // input reader calls FillComplete(). If you create your own discretizations
  // try to understand what you do.

  // Get highest GID used so far and add one
  int count = GetFirstGIDNumberToBeUsed(dis);

  // Check if we have a face discretization which supports degrees of freedom on faces
  Teuchos::RCP<const DiscretizationHDG> facedis =
      Teuchos::rcp_dynamic_cast<const DiscretizationHDG>(Teuchos::rcp(&dis, false));

  // set count to 0 in case of dofset 2 in HDG discretizations
  if (facedis != Teuchos::null && dspos_ == 2) count = 0;

  // Now this is tricky. We have to care for nodes, faces, and elements, both
  // row and column maps. In general both nodes, faces, and elements can have
  // dofs. In all cases these dofs might be shared with other nodes, faces,
  // or elements. (The very general case. For elements we'd probably
  // don't need that.)
  //
  // The point is that we have to make sure the dof numbering of a
  // mesh is independent of its parallel distribution. Otherwise we
  // could not redistribute a mesh. We would not be able to use old
  // distributed vectors afterwards.
  //
  // Each object (node or element) could have a different number of
  // dofs. The parallel distribution is arbitrary. So we fall back to
  // two redundant vectors here to gather the number of dofs per node
  // or element.

  // numdf for all nodes and elements
  numdfcolnodes_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
  numdfcolelements_ = Teuchos::rcp(new Epetra_IntVector(*dis.ElementColMap()));
  if (facedis != Teuchos::null && facedis->FaceColMap() != nullptr)
    numdfcolfaces_ = Teuchos::rcp(new Epetra_IntVector(*facedis->FaceColMap()));

  // index of first dof for all nodes and elements
  idxcolnodes_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
  idxcolelements_ = Teuchos::rcp(new Epetra_IntVector(*dis.ElementColMap()));
  if (facedis != Teuchos::null && facedis->FaceColMap() != nullptr)
    idxcolfaces_ = Teuchos::rcp(new Epetra_IntVector(*facedis->FaceColMap()));

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  int maxnodenumdf = 0;
  int maxelementnumdf = 0;
  std::map<int, std::vector<int>> nodedofset;
  std::map<int, std::vector<int>> nodeduplicatedofset;
  std::map<int, std::vector<int>> elementdofset;
  std::map<int, std::vector<int>> facedofset;

  {
    // get DoF coupling conditions
    std::vector<DRT::Condition*> couplingconditions(0);
    dis.GetCondition("PointCoupling", couplingconditions);
    if ((int)couplingconditions.size() > 0) pccdofhandling_ = true;

    // do the nodes first
    Epetra_IntVector numdfrownodes(*dis.NodeRowMap());
    Epetra_IntVector idxrownodes(*dis.NodeRowMap());

    int numrownodes = dis.NumMyRowNodes();
    for (int i = 0; i < numrownodes; ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      numdfrownodes[i] = NumDofPerNode(*actnode);
    }

    int minnodegid = GetMinimalNodeGIDIfRelevant(dis);
    maxnodenumdf = numdfrownodes.MaxValue();
    GetReservedMaxNumDofperNode(maxnodenumdf);  // XFEM::XFEMDofSet set to const number!

    for (int i = 0; i < numrownodes; ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      const int gid = actnode->Id();

      // **********************************************************************
      // **********************************************************************
      // check for DoF coupling conditions                         popp 02/2016
      // **********************************************************************
      // **********************************************************************
      int relevantcondid = -1;
      if (dspos_ == 0)
      {
        for (int k = 0; k < (int)couplingconditions.size(); ++k)
        {
          if (couplingconditions[k]->ContainsNode(gid))
          {
            if (relevantcondid != -1) dserror("ERROR: Two coupling conditions on one node");
            relevantcondid = k;
          }
        }
      }

      // check for node coupling condition and slave/master status
      bool specialtreatment = false;
      if (relevantcondid >= 0)
      {
        const std::vector<int>* nodeids = couplingconditions[relevantcondid]->Nodes();
        if (!nodeids) dserror("ERROR: Condition does not have Node Ids");

        // check if all nodes in this condition are on same processor
        // (otherwise throw a dserror for now - not yet implemented)
        bool allononeproc = true;
        for (int k = 0; k < (int)(nodeids->size()); ++k)
        {
          int checkgid = (*nodeids)[k];
          if (!dis.NodeRowMap()->MyGID(checkgid)) allononeproc = false;
        }
        if (!allononeproc)
          dserror(
              "ERRROR: Nodes in point coupling condition must all be on same processsor (for now)");

        // do nothing for first (master) node in coupling condition
        // do something for second, third, ... (slave) node
        if ((*nodeids)[0] != gid)
        {
          // critical case
          specialtreatment = true;

          // check total number of dofs and determine which dofs are to be coupled
          if (couplingconditions[relevantcondid]->GetInt("numdof") != numdfrownodes[i])
            dserror("ERROR: Number of DoFs in coupling condition (%i) does not match node (%i)",
                couplingconditions[relevantcondid]->GetInt("numdof"), numdfrownodes[i]);
          const std::vector<int>* onoffcond =
              couplingconditions[relevantcondid]->Get<std::vector<int>>("onoff");

          // get master node of this condition
          int mgid = (*nodeids)[0];
          std::vector<int>& mdofs = nodedofset[mgid];
          if ((int)(mdofs.size()) == 0)
            dserror("ERROR: Master node has not yet been initialized with DoFs");

          // special treatment
          int numdf = numdfrownodes[i];
          int dof = count + (gid - minnodegid) * maxnodenumdf;
          idxrownodes[i] = dof;
          std::vector<int>& dofs = nodedofset[gid];
          std::vector<int>& duplicatedofs = nodeduplicatedofset[gid];
          dofs.reserve(numdf);
          duplicatedofs.reserve(numdf);
          for (int j = 0; j < numdf; ++j)
          {
            // push back master node DoF ID if coupled
            if ((*onoffcond)[j] == 1)
            {
              dofs.push_back(mdofs[j]);
              duplicatedofs.push_back(1);
            }
            // push back new DoF ID if not coupled
            else
            {
              dofs.push_back(dof + j);
              duplicatedofs.push_back(0);
            }
          }
        }
      }

      // standard treatment for non-coupling nodes and master coupling nodes
      if (!specialtreatment)
      {
        int numdf = numdfrownodes[i];
        int dof = count + (gid - minnodegid) * maxnodenumdf;
        idxrownodes[i] = dof;
        std::vector<int>& dofs = nodedofset[gid];
        std::vector<int>& duplicatedofs = nodeduplicatedofset[gid];
        dofs.reserve(numdf);
        duplicatedofs.reserve(numdf);
        for (int j = 0; j < numdf; ++j)
        {
          dofs.push_back(dof + j);
          duplicatedofs.push_back(0);
        }
      }
      // **********************************************************************
      // **********************************************************************
      // **********************************************************************
      // **********************************************************************
    }

    Epetra_Import nodeimporter(numdfcolnodes_->Map(), numdfrownodes.Map());
    int err = numdfcolnodes_->Import(numdfrownodes, nodeimporter, Insert);
    if (err) dserror("Import using importer returned err=%d", err);
    err = idxcolnodes_->Import(idxrownodes, nodeimporter, Insert);
    if (err) dserror("Import using importer returned err=%d", err);

    count = maxnodenumdf > 0 ? idxrownodes.MaxValue() + maxnodenumdf : 0;

    //////////////////////////////////////////////////////////////////

    // Now do it again for the faces
    if (facedis != Teuchos::null && facedis->FaceRowMap() != nullptr)
    {
      Epetra_IntVector numdfrowfaces(*facedis->FaceRowMap());
      Epetra_IntVector idxrowfaces(*facedis->FaceRowMap());
      int numcolelements = dis.NumMyColElements();

      const int mypid = dis.Comm().MyPID();
      for (int i = 0; i < numcolelements; ++i)
      {
        Teuchos::RCP<DRT::FaceElement>* faces = dis.lColElement(i)->Faces();
        // If no faces are found, continue...
        if (faces == nullptr) continue;
        for (int face = 0; face < dis.lColElement(i)->NumFace(); ++face)
          if (faces[face]->Owner() == mypid)
          {
            const int mylid = facedis->FaceRowMap()->LID(faces[face]->Id());
            numdfrowfaces[mylid] = NumDofPerFace(*(dis.lColElement(i)), face);
          }
      }

      int minfacegid = facedis->FaceRowMap()->MinAllGID();
      int maxfacenumdf = numdfrowfaces.MaxValue();

      for (int i = 0; i < numcolelements; ++i)
      {
        Teuchos::RCP<DRT::FaceElement>* faces = dis.lColElement(i)->Faces();
        if (faces == nullptr) continue;
        for (int face = 0; face < dis.lColElement(i)->NumFace(); ++face)
          if (faces[face]->Owner() == mypid)
          {
            const int gid = faces[face]->Id();
            const int mylid = facedis->FaceRowMap()->LID(gid);
            int numdf = numdfrowfaces[mylid];
            int dof = count + (gid - minfacegid) * maxfacenumdf;
            idxrowfaces[mylid] = dof;
            std::vector<int>& dofs = facedofset[gid];
            // do not visit the same face more than once
            if (dofs.empty())
            {
              dofs.reserve(numdf);
              for (int j = 0; j < numdf; ++j)
              {
                dofs.push_back(dof + j);
              }
            }
          }
      }

      Epetra_Import faceimporter(numdfcolfaces_->Map(), numdfrowfaces.Map());
      err = numdfcolfaces_->Import(numdfrowfaces, faceimporter, Insert);
      if (err) dserror("Import using importer returned err=%d", err);
      err = idxcolfaces_->Import(idxrowfaces, faceimporter, Insert);
      if (err) dserror("Import using importer returned err=%d", err);

      count = idxrowfaces.MaxValue() + maxfacenumdf;
    }

    //////////////////////////////////////////////////////////////////

    // Now do it again for the elements
    Epetra_IntVector numdfrowelements(*dis.ElementRowMap());
    Epetra_IntVector idxrowelements(*dis.ElementRowMap());

    int numrowelements = dis.NumMyRowElements();
    for (int i = 0; i < numrowelements; ++i)
    {
      DRT::Element* actele = dis.lRowElement(i);
      // const int gid = actele->Id();
      int numdf = NumDofPerElement(*actele);
      numdfrowelements[i] = numdf;
    }

    int minelementgid = dis.ElementRowMap()->MinAllGID();
    maxelementnumdf = numdfrowelements.MaxValue();

    for (int i = 0; i < numrowelements; ++i)
    {
      DRT::Element* actelement = dis.lRowElement(i);
      const int gid = actelement->Id();
      int numdf = numdfrowelements[i];
      int dof = count + (gid - minelementgid) * maxelementnumdf;
      idxrowelements[i] = dof;
      std::vector<int>& dofs = elementdofset[gid];
      dofs.reserve(numdf);
      for (int j = 0; j < numdf; ++j)
      {
        dofs.push_back(dof + j);
      }
    }

    Epetra_Import elementimporter(numdfcolelements_->Map(), numdfrowelements.Map());
    err = numdfcolelements_->Import(numdfrowelements, elementimporter, Insert);
    if (err) dserror("Import using importer returned err=%d", err);
    err = idxcolelements_->Import(idxrowelements, elementimporter, Insert);
    if (err) dserror("Import using importer returned err=%d", err);
  }

  // Now finally we have everything in place to build the maps.
  int numrownodes = dis.NumMyRowNodes();
  int numrowelements = dis.NumMyRowElements();

  std::vector<int> localrowdofs;
  std::vector<int> localcoldofs;
  localrowdofs.reserve(numrownodes * maxnodenumdf + numrowelements * maxelementnumdf);
  localcoldofs.reserve(numrownodes * maxnodenumdf + numrowelements * maxelementnumdf);

  std::vector<int> allnodelocalcoldofs;
  allnodelocalcoldofs.reserve(numrownodes * maxnodenumdf);

  for (std::map<int, std::vector<int>>::iterator i = nodedofset.begin(); i != nodedofset.end(); ++i)
  {
    std::vector<int>& dofs = i->second;
    std::vector<int>& duplicatedofs = nodeduplicatedofset[i->first];
    std::vector<int> cleandofs;
    for (unsigned j = 0; j < dofs.size(); ++j)
    {
      if (duplicatedofs[j] == 0) cleandofs.push_back(dofs[j]);
    }
    std::copy(cleandofs.begin(), cleandofs.end(), std::back_inserter(localrowdofs));
    // printf("Proc %d nodal gid %d ndofs %d\n",proc,i->first,(int)dofs.size());
    // for (unsigned j=0; j<dofs.size(); ++j) printf(" %d ",dofs[j]);
    // printf("\n");
  }
  for (std::map<int, std::vector<int>>::iterator i = facedofset.begin(); i != facedofset.end(); ++i)
  {
    std::vector<int>& dofs = i->second;
    std::copy(dofs.begin(), dofs.end(), std::back_inserter(localrowdofs));
    // printf("Proc %d ele gid %d ndofs %d\n",dis.Comm().MyPID(),i->first,(int)dofs.size());
    // for (unsigned j=0; j<dofs.size(); ++j) printf(" %d ",dofs[j]);
    // printf("\n");
  }
  for (std::map<int, std::vector<int>>::iterator i = elementdofset.begin();
       i != elementdofset.end(); ++i)
  {
    std::vector<int>& dofs = i->second;
    std::copy(dofs.begin(), dofs.end(), std::back_inserter(localrowdofs));
    // printf("Proc %d ele gid %d ndofs %d\n",dis.Comm().MyPID(),i->first,(int)dofs.size());
    // for (unsigned j=0; j<dofs.size(); ++j) printf(" %d ",dofs[j]);
    // printf("\n");
  }

  CORE::COMM::Exporter nodeexporter(*dis.NodeRowMap(), *dis.NodeColMap(), dis.Comm());
  nodeexporter.Export(nodedofset);

  CORE::COMM::Exporter elementexporter(*dis.ElementRowMap(), *dis.ElementColMap(), dis.Comm());
  elementexporter.Export(elementdofset);

  if (facedis != Teuchos::null && facedis->FaceRowMap() != nullptr)
  {
    CORE::COMM::Exporter faceexporter(*facedis->FaceRowMap(), *facedis->FaceColMap(), dis.Comm());
    faceexporter.Export(facedofset);
  }

  for (std::map<int, std::vector<int>>::iterator i = nodedofset.begin(); i != nodedofset.end(); ++i)
  {
    std::vector<int>& dofs = i->second;
    std::vector<int> cleandofs;
    for (unsigned j = 0; j < dofs.size(); ++j)
    {
      if (std::find(localcoldofs.begin(), localcoldofs.end(), dofs[j]) == localcoldofs.end())
        cleandofs.push_back(dofs[j]);
    }
    std::copy(cleandofs.begin(), cleandofs.end(), std::back_inserter(localcoldofs));
    std::copy(dofs.begin(), dofs.end(), std::back_inserter(allnodelocalcoldofs));
  }
  for (std::map<int, std::vector<int>>::iterator i = facedofset.begin(); i != facedofset.end(); ++i)
  {
    std::vector<int>& dofs = i->second;
    std::copy(dofs.begin(), dofs.end(), std::back_inserter(localcoldofs));
  }
  for (std::map<int, std::vector<int>>::iterator i = elementdofset.begin();
       i != elementdofset.end(); ++i)
  {
    std::vector<int>& dofs = i->second;
    std::copy(dofs.begin(), dofs.end(), std::back_inserter(localcoldofs));
  }

  dofrowmap_ =
      Teuchos::rcp(new Epetra_Map(-1, localrowdofs.size(), localrowdofs.data(), 0, dis.Comm()));
  if (!dofrowmap_->UniqueGIDs()) dserror("Dof row map is not unique");
  dofcolmap_ =
      Teuchos::rcp(new Epetra_Map(-1, localcoldofs.size(), localcoldofs.data(), 0, dis.Comm()));

  // **********************************************************************
  // **********************************************************************
  // build map of all (non-unique) column DoFs
  dofscolnodes_ = Teuchos::rcp(
      new Epetra_Map(-1, allnodelocalcoldofs.size(), allnodelocalcoldofs.data(), 0, dis.Comm()));

  // build shift vector
  shiftcolnodes_ = Teuchos::rcp(new Epetra_IntVector(*dis.NodeColMap()));
  int numcolnodes = dis.NumMyColNodes();
  for (int i = 0; i < numcolnodes; ++i)
  {
    if (i == 0)
    {
      (*shiftcolnodes_)[i] = 0;
    }
    else
    {
      DRT::Node* lastnode = dis.lColNode(i - 1);
      (*shiftcolnodes_)[i] = (*shiftcolnodes_)[i - 1] + NumDofPerNode(*lastnode);
    }
  }
  // **********************************************************************
  // **********************************************************************

  // degrees of freedom have now been assigned
  filled_ = true;

  // tell all proxies
  NotifyAssigned();

  // return maximum dof number of this dofset (+1)
  count = dofrowmap_->MaxAllGID() + 1;
  return count;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::DofSet::Initialized() const
{
  if (dofcolmap_ == Teuchos::null or dofrowmap_ == Teuchos::null)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DofSet::DofRowMap() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSet::DofRowMap(): dofrowmap_ not initialized, yet");
  return dofrowmap_.get();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DofSet::DofColMap() const
{
  if (dofcolmap_ == Teuchos::null)
    dserror("DRT::DofSet::DofColMap(): dofcolmap_ not initialized, yet");
  return dofcolmap_.get();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSet::NumGlobalElements() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSet::NumGlobalElements(): dofrowmap_ not initialized, yet");
  return dofrowmap_->NumGlobalElements();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSet::MaxAllGID() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSet::MaxAllGID(): dofrowmap_ not initialized, yet");
  return dofrowmap_->MaxAllGID();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSet::MinAllGID() const
{
  if (dofrowmap_ == Teuchos::null)
    dserror("DRT::DofSet::MinAllGID(): dofrowmap_ not initialized, yet");
  return dofrowmap_->MinAllGID();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSet::GetFirstGIDNumberToBeUsed(const Discretization& dis) const
{
  return MaxGIDinList(dis.Comm()) + 1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSet::GetMinimalNodeGIDIfRelevant(const Discretization& dis) const
{
  return dis.NodeRowMap()->MinAllGID();
}

BACI_NAMESPACE_CLOSE
