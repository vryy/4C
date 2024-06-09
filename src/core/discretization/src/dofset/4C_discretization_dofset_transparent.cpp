/*---------------------------------------------------------------------*/
/*! \file

\brief A set of degrees of freedom special for contact

\level 1


*/
/*---------------------------------------------------------------------*/

#include "4C_discretization_dofset_transparent.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN

Core::DOFSets::TransparentDofSet::TransparentDofSet(
    Teuchos::RCP<Discret::Discretization> sourcedis, bool parallel)
    : Core::DOFSets::DofSet(), sourcedis_(sourcedis), parallel_(parallel)
{
  return;
}

int Core::DOFSets::TransparentDofSet::assign_degrees_of_freedom(
    const Discret::Discretization& dis, const unsigned dspos, const int start)
{
  // first, we call the standard assign_degrees_of_freedom from the base class
  int count = DofSet::assign_degrees_of_freedom(dis, dspos, start);
  if (pccdofhandling_)
    FOUR_C_THROW("ERROR: Point coupling cinditions not yet implemented for TransparentDofSet");

  if (!parallel_)
  {
    transfer_degrees_of_freedom(*sourcedis_, dis, start);
  }
  else
  {
    parallel_transfer_degrees_of_freedom(*sourcedis_, dis, start);
  }

  // tell all proxies (again!)
  NotifyAssigned();

  return count;
}

/// Assign dof numbers for new discretization using dof numbering from source discretization.
void Core::DOFSets::TransparentDofSet::transfer_degrees_of_freedom(
    const Discret::Discretization& sourcedis, const Discret::Discretization& newdis,
    const int start)
{
  if (!sourcedis.dof_row_map()->UniqueGIDs()) FOUR_C_THROW("dof_row_map is not unique");
  if (!sourcedis.NodeRowMap()->UniqueGIDs()) FOUR_C_THROW("NodeRowMap is not unique");
  if (!sourcedis.ElementRowMap()->UniqueGIDs()) FOUR_C_THROW("ElementRowMap is not unique");

  if (!newdis.dof_row_map()->UniqueGIDs()) FOUR_C_THROW("dof_row_map is not unique");
  if (!newdis.NodeRowMap()->UniqueGIDs()) FOUR_C_THROW("NodeRowMap is not unique");
  if (!newdis.ElementRowMap()->UniqueGIDs()) FOUR_C_THROW("ElementRowMap is not unique");

  // build dofrowmap
  std::set<int> dofrowset;
  std::vector<int> dofrowvec;
  dofrowvec.reserve(dofrowmap_->NumMyElements());
  for (int inode = 0; inode != newdis.NumMyRowNodes(); ++inode)
  {
    const Core::Nodes::Node* newnode = newdis.lRowNode(inode);
    const Core::Nodes::Node* sourcenode = sourcedis.gNode(newnode->Id());

    const std::vector<int> dofs = sourcedis.Dof(0, sourcenode);

    const int newlid = newnode->LID();
    const int numdofs = (*numdfcolnodes_)[newlid];
    if (numdofs > 0)
    {
      (*idxcolnodes_)[newlid] = dofs[0];
      for (int idof = 0; idof < numdofs; ++idof)
      {
        dofrowset.insert(dofs[idof]);
      }
    }
  }

  for (std::set<int>::iterator idof = dofrowset.begin(); idof != dofrowset.end(); ++idof)
  {
    dofrowvec.push_back(*idof);
  }

  dofrowmap_ =
      Teuchos::rcp(new Epetra_Map(-1, dofrowvec.size(), dofrowvec.data(), 0, newdis.Comm()));

  // build dofcolvec
  std::set<int> dofcolset;
  std::vector<int> dofcolvec;
  dofcolvec.reserve(dofcolmap_->NumMyElements());
  for (int inode = 0; inode != newdis.NumMyColNodes(); ++inode)
  {
    const Core::Nodes::Node* newnode = newdis.lColNode(inode);
    const Core::Nodes::Node* sourcenode = sourcedis.gNode(newnode->Id());

    const int lid = sourcenode->LID();
    if (lid == -1)
    {
      FOUR_C_THROW("required node %d not on proc", newnode->Id());
    }
    const std::vector<int> dofs = sourcedis.Dof(0, sourcenode);
    const int newlid = newnode->LID();
    // const int newfirstidx = (*idxcolnodes_)[newlid];
    const int numdofs = (*numdfcolnodes_)[newlid];
    if (numdofs > 0)
    {
      (*idxcolnodes_)[newlid] = dofs[0];
      //        if(numdofs!=(int)dofs.size())
      //        FOUR_C_THROW("numdofs %d!=%d for node %d",numdofs,(int)dofs.size(),newnode->Id());

      for (int idof = 0; idof < numdofs; ++idof)
      {
        dofcolset.insert(dofs[idof]);
      }
    }
  }

  for (std::set<int>::iterator idof = dofcolset.begin(); idof != dofcolset.end(); ++idof)
  {
    dofcolvec.push_back(*idof);
  }

  dofcolmap_ =
      Teuchos::rcp(new Epetra_Map(-1, dofcolvec.size(), dofcolvec.data(), 0, newdis.Comm()));
}

/// Assign dof numbers for new discretization using dof numbering from source discretization.
void Core::DOFSets::TransparentDofSet::parallel_transfer_degrees_of_freedom(
    const Discret::Discretization& sourcedis, const Discret::Discretization& newdis,
    const int start)
{
  if (!sourcedis.dof_row_map()->UniqueGIDs()) FOUR_C_THROW("dof_row_map is not unique");
  if (!sourcedis.NodeRowMap()->UniqueGIDs()) FOUR_C_THROW("NodeRowMap is not unique");
  if (!sourcedis.ElementRowMap()->UniqueGIDs()) FOUR_C_THROW("ElementRowMap is not unique");

  if (!newdis.dof_row_map()->UniqueGIDs()) FOUR_C_THROW("dof_row_map is not unique");
  if (!newdis.NodeRowMap()->UniqueGIDs()) FOUR_C_THROW("NodeRowMap is not unique");
  if (!newdis.ElementRowMap()->UniqueGIDs()) FOUR_C_THROW("ElementRowMap is not unique");

  // list all my rownode ids

  //
  // we need a mapping
  //
  // colnode gid -> std::vector<int> dofs of sourcenode
  //
  // problem: sourcenode not necessarily on this proc -> communicate
  //
  // the idea is to search for the sourcerownode on some proc and to get
  // this unique number
  //
  std::map<int, std::vector<int>> gid_to_dofs;

  for (int inode = 0; inode != newdis.NumMyColNodes(); ++inode)
  {
    const Core::Nodes::Node* newnode = newdis.lColNode(inode);
    int gid = newnode->Id();
    std::vector<int> emptyvec;
    gid_to_dofs.insert(std::pair<int, std::vector<int>>(gid, emptyvec));
  }

  {
    // create an exporter for point to point comunication
    Core::Communication::Exporter exporter(sourcedis.Comm());

    // necessary variables
    MPI_Request request;

    // define send and receive blocks
    std::vector<char> sblock;
    std::vector<char> rblock;

    // get number of processors and the current processors id
    int numproc = sourcedis.Comm().NumProc();
    int myrank = sourcedis.Comm().MyPID();

    //----------------------------------------------------------------------
    // communication is done in a round robin loop
    //----------------------------------------------------------------------
    for (int np = 0; np < numproc + 1; ++np)
    {
      // in the first step, we cannot receive anything
      if (np > 0)
      {
        receive_block(numproc, myrank, rblock, exporter, request);

        // Unpack info from the receive block from the last proc
        unpack_local_source_dofs(gid_to_dofs, rblock);
      }

      // in the last step, we keep everything on this proc
      if (np < numproc)
      {
        // -----------------------
        // do what we wanted to do
        set_source_dofs_available_on_this_proc(gid_to_dofs);

        // Pack info into block to send
        Core::Communication::PackBuffer data;
        PackLocalSourceDofs(gid_to_dofs, data);
        gid_to_dofs.clear();
        swap(sblock, data());

        send_block(numproc, myrank, sblock, exporter, request);
      }
    }
  }

  std::set<int> slaveset;
  std::vector<Core::Conditions::Condition*> mypbcs;

  // get periodic surface boundary conditions
  sourcedis_->GetCondition("SurfacePeriodic", mypbcs);

  if (mypbcs.empty())
  {
    sourcedis_->GetCondition("LinePeriodic", mypbcs);
  }

  for (unsigned numcond = 0; numcond < mypbcs.size(); ++numcond)
  {
    Core::Conditions::Condition* thiscond = mypbcs[numcond];

    // see whether we have a slave condition
    const std::string& mymasterslavetoggle =
        thiscond->parameters().Get<std::string>("Is slave periodic boundary condition");

    if (!(mymasterslavetoggle == "Master"))
    {
      const std::vector<int>* pbcids;
      pbcids = (*thiscond).GetNodes();

      for (std::vector<int>::const_iterator iter = pbcids->begin(); iter != pbcids->end(); ++iter)
      {
        slaveset.insert(*iter);
      }
    }
  }

  // build dofrowmap
  std::set<int> dofrowset;
  std::vector<int> dofrowvec;
  dofrowvec.reserve(dofrowmap_->NumMyElements());
  for (int inode = 0; inode != newdis.NumMyRowNodes(); ++inode)
  {
    const Core::Nodes::Node* newnode = newdis.lRowNode(inode);

    const std::vector<int> dofs = gid_to_dofs[newnode->Id()];

    const int newlid = newnode->LID();
    const int numdofs = (*numdfcolnodes_)[newlid];

    if (numdofs != (int)dofs.size())
    {
      printf("This is node %d  (%12.5e,%12.5e,%12.5e)\n", newnode->Id(), newnode->X()[0],
          newnode->X()[1], newnode->X()[2]);

      FOUR_C_THROW("spooky, isn't it? dofs to overwrite %d != %d dofs.size() to set \n", numdofs,
          dofs.size());
    }

    if (numdofs > 0)
    {
      (*idxcolnodes_)[newlid] = dofs[0];

      // slave-dofs must not enter the dofrowset (if master&slave are on different procs)
      std::set<int>::iterator curr = slaveset.find(newnode->Id());

      if (curr == slaveset.end())
      {
        for (int idof = 0; idof < numdofs; ++idof)
        {
          dofrowset.insert(dofs[idof]);
        }
      }
    }
  }

  for (std::set<int>::iterator idof = dofrowset.begin(); idof != dofrowset.end(); ++idof)
  {
    dofrowvec.push_back(*idof);
  }

  dofrowmap_ =
      Teuchos::rcp(new Epetra_Map(-1, dofrowvec.size(), dofrowvec.data(), 0, newdis.Comm()));

  // build dofcolvec
  std::set<int> dofcolset;
  std::vector<int> dofcolvec;
  dofcolvec.reserve(dofcolmap_->NumMyElements());
  for (int inode = 0; inode != newdis.NumMyColNodes(); ++inode)
  {
    const Core::Nodes::Node* newnode = newdis.lColNode(inode);

    const std::vector<int> dofs = gid_to_dofs[newnode->Id()];

    const int newlid = newnode->LID();
    // const int newfirstidx = (*idxcolnodes_)[newlid];
    const int numdofs = (*numdfcolnodes_)[newlid];
    if (numdofs > 0)
    {
      (*idxcolnodes_)[newlid] = dofs[0];
      if (numdofs != (int)dofs.size())
        FOUR_C_THROW("numdofs %d!=%d for node %d", numdofs, dofs.size(), newnode->Id());

      for (int idof = 0; idof < numdofs; ++idof)
      {
        dofcolset.insert(dofs[idof]);
      }
    }
  }

  for (std::set<int>::iterator idof = dofcolset.begin(); idof != dofcolset.end(); ++idof)
  {
    dofcolvec.push_back(*idof);
  }

  dofcolmap_ =
      Teuchos::rcp(new Epetra_Map(-1, dofcolvec.size(), dofcolvec.data(), 0, newdis.Comm()));


  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | write dof info into map                                    (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Core::DOFSets::TransparentDofSet::set_source_dofs_available_on_this_proc(
    std::map<int, std::vector<int>>& gid_to_dofs)
{
  for (std::map<int, std::vector<int>>::iterator curr = gid_to_dofs.begin();
       curr != gid_to_dofs.end(); ++curr)
  {
    const int lid = sourcedis_->NodeRowMap()->LID(curr->first);

    if (lid > -1)
    {
      curr->second.clear();

      const Core::Nodes::Node* sourcenode = sourcedis_->gNode(curr->first);

      const std::vector<int> dofs = sourcedis_->Dof(0, sourcenode);

      for (std::vector<int>::const_iterator iter = dofs.begin(); iter != dofs.end(); ++iter)
      {
        curr->second.push_back(*iter);
      }
    }
    else
    {
      int numproc = sourcedis_->Comm().NumProc();
      if (numproc == 1)
      {
        FOUR_C_THROW(
            "I have a one-processor problem but the node is not on the proc. "
            "sourcedis_->NodeRowMap() is probably currupted.");
      }
    }
  }
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | pack all values into a send block                          (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Core::DOFSets::TransparentDofSet::PackLocalSourceDofs(
    std::map<int, std::vector<int>>& gid_to_dofs, Core::Communication::PackBuffer& sblock)
{
  int size = gid_to_dofs.size();

  // add size  to sendblock
  Core::Communication::ParObject::add_to_pack(sblock, size);

  for (std::map<int, std::vector<int>>::iterator curr = gid_to_dofs.begin();
       curr != gid_to_dofs.end(); ++curr)
  {
    int gid = curr->first;
    std::vector<int> mydofs = curr->second;
    int numdofs = (int)mydofs.size();

    Core::Communication::ParObject::add_to_pack(sblock, gid);
    Core::Communication::ParObject::add_to_pack(sblock, numdofs);
    for (int ll = 0; ll < numdofs; ++ll)
    {
      Core::Communication::ParObject::add_to_pack(sblock, mydofs[ll]);
    }
  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | unpack all values contained in receive block               (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Core::DOFSets::TransparentDofSet::unpack_local_source_dofs(
    std::map<int, std::vector<int>>& gid_to_dofs, std::vector<char>& rblock)
{
  gid_to_dofs.clear();

  // position to extract
  std::vector<char>::size_type position = 0;

  // extract size
  int size = 0;
  Core::Communication::ParObject::extract_from_pack(position, rblock, size);

  for (int rr = 0; rr < size; ++rr)
  {
    int gid = -1;
    std::vector<int> mydofs;
    int numdofs = 0;

    Core::Communication::ParObject::extract_from_pack(position, rblock, gid);
    Core::Communication::ParObject::extract_from_pack(position, rblock, numdofs);

    for (int ll = 0; ll < numdofs; ++ll)
    {
      int thisdof = 0;

      Core::Communication::ParObject::extract_from_pack(position, rblock, thisdof);
      mydofs.push_back(thisdof);
    }

    gid_to_dofs.insert(std::pair<int, std::vector<int>>(gid, mydofs));
  }

  rblock.clear();
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | receive a block in the round robin communication pattern   (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Core::DOFSets::TransparentDofSet::receive_block(int numproc, int myrank,
    std::vector<char>& rblock, Core::Communication::Exporter& exporter, MPI_Request& request)
{
  // necessary variables

  int length = -1;
  int frompid = (myrank + numproc - 1) % numproc;
  int tag = frompid;

  // make sure that you do not think you received something if
  // you didn't
  if (rblock.empty() == false)
  {
    FOUR_C_THROW("rblock not empty");
  }

  // receive from predecessor
  exporter.ReceiveAny(frompid, tag, rblock, length);

  if (tag != (myrank + numproc - 1) % numproc)
  {
    FOUR_C_THROW("received wrong message (ReceiveAny)");
  }

  exporter.Wait(request);

  // for safety
  exporter.Comm().Barrier();

  return;
}  // receive_block


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | send a block in the round robin communication pattern      (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Core::DOFSets::TransparentDofSet::send_block(int numproc, int myrank,
    std::vector<char>& sblock, Core::Communication::Exporter& exporter, MPI_Request& request)
{
  // Send block to next proc.
  int tag = myrank;
  int frompid = myrank;
  int topid = (myrank + 1) % numproc;

  exporter.i_send(frompid, topid, sblock.data(), sblock.size(), tag, request);


  // for safety
  exporter.Comm().Barrier();

  return;
}  // send_block

FOUR_C_NAMESPACE_CLOSE
