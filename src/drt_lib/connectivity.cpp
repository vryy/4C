/*----------------------------------------------------------------------------*/
/*! \file

\brief Find connectivities in a discretization.

\maintainer Martin Kronbichler

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "connectivity.H"

#include "../drt_lib/drt_discret_interface.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_exporter.H"

#include <Epetra_Map.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::Connectivity::splitDiscretizationIntoDistinctNodeMaps(
    const DRT::DiscretizationInterface& full_discret, const MapFormat format,
    std::vector<Teuchos::RCP<Epetra_Map>>& body_node_maps)
{
  body_node_maps.clear();

  node_sets my_body_nsets;
  std::vector<int> my_new_body_ids;

  // interface nodes
  const int num_inodes = full_discret.NodeRowMap()->NumMyElements();  // inmap.NumMyElements();
  const int* ingids = full_discret.NodeRowMap()->MyGlobalElements();

  for (int i = 0; i < num_inodes; ++i)
  {
    // Take one node from the interface discretization and find it in the
    // full discretization.
    const int ingid = ingids[i];
    DRT::Node* node = full_discret.gNode(ingid);

    // test if this node is already in one of our node sets
    unsigned body_id = 0;
    for (auto& my_body_nset : my_body_nsets)
    {
      if (my_body_nset.find(node) != my_body_nset.end()) break;
      ++body_id;
    }

    // if we could find it, there is nothing more to do since it has already
    // been considered
    if (body_id < my_body_nsets.size()) continue;

    // if we can't find the current node in one of the already
    // processed node sets, it must belong to a new body/set
    my_new_body_ids.push_back(my_body_nsets.size());
    my_body_nsets.push_back(node_set());
    node_set& my_body_nset = my_body_nsets.back();
    my_body_nset.insert(node);
    fillMyConnectedNodeSet(my_body_nset, full_discret);
  }

  const Epetra_Comm& comm = full_discret.Comm();

  int my_num_new_bodies = my_new_body_ids.size();
  int g_num_new_bodies = 0;
  comm.SumAll(&my_num_new_bodies, &g_num_new_bodies, 1);

  Connection connections(comm.MyPID());

  while (g_num_new_bodies > 0)
  {
    int my_num_bodies = my_body_nsets.size();
    int max_num_bodies = my_body_nsets.size();
    comm.MaxAll(&my_num_bodies, &max_num_bodies, 1);

    std::vector<node_set*> my_new_body_nsets(max_num_bodies, NULL);
    for (int my_new_body_id : my_new_body_ids)
      my_new_body_nsets[my_new_body_id] = &my_body_nsets[my_new_body_id];

    received_gid_set all_received_ngids(max_num_bodies);
    receiveAnyOfMyGids(all_received_ngids, my_new_body_nsets, comm);

    // body IDs which are connected to already existing bodies on other procs,
    // but which were not found on this proc in the previous attempt
    my_new_body_ids.clear();

    // loop over all received bodies
    for (unsigned rec_b = 0; rec_b < all_received_ngids.size(); ++rec_b)
    {
      auto& received_ngids = all_received_ngids[rec_b];

      // received nodal gids from the KEY proc
      for (auto& ngids_fp : received_ngids)
      {
        const int from_proc = ngids_fp.first;
        // looped over the ghosted gids of proc "from_proc" and find the correct
        // body on this proc which is the owner of the corresponding nodes
        for (int mygid : ngids_fp.second)
        {
          // Take one node from the interface discretization and find it in the
          // full discretization.
          DRT::Node* node = full_discret.gNode(mygid);

          // test if this node is already in one of our node sets
          unsigned body_id = 0;
          for (auto& my_body_nset : my_body_nsets)
          {
            if (my_body_nset.find(node) != my_body_nset.end()) break;
            ++body_id;
          }

          // could find a body
          if (body_id < my_body_nsets.size())
          {
            connections.add(body_id, rec_b, from_proc);
            break;
          }

          // since we couldn't find the current node in one of the already
          // processed node sets, it must belong to a new body/set on this proc
          my_new_body_ids.push_back(my_body_nsets.size());
          my_body_nsets.push_back(node_set());
          node_set& my_body_nset = my_body_nsets.back();
          my_body_nset.insert(node);
          fillMyConnectedNodeSet(my_body_nset, full_discret);
        }
      }
      //      for ( int new_id : my_new_body_ids )
      //        std::cout << "New body with id #" << new_id << " found on proc #"
      //            << comm.MyPID() << std::endl;
    }
    my_num_new_bodies = my_new_body_ids.size();
    g_num_new_bodies = 0;
    comm.SumAll(&my_num_new_bodies, &g_num_new_bodies, 1);
  }

  //  connections.print(std::cout);
  createBodyNodeMaps(format, comm, connections, my_body_nsets, body_node_maps);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::Connectivity::createBodyNodeMaps(const MapFormat format, const Epetra_Comm& comm,
    const Connection& connections, const node_sets& my_body_nsets,
    std::vector<Teuchos::RCP<Epetra_Map>>& body_node_maps)
{
  DRT::Exporter exporter(comm);

  const int num_procs = comm.NumProc();

  std::vector<int> sendto(num_procs, 0);
  std::vector<int> receivefrom(num_procs, 0);
  std::vector<MPI_Request> mpi_requests;

  std::set<int> done;
  int p = 0;
  while (p < num_procs)
  {
    int build_bid = -1;
    int next_proc = p;
    std::fill(sendto.data(), sendto.data() + num_procs, 0);
    std::fill(receivefrom.data(), receivefrom.data() + num_procs, 0);

    std::set<int> already_considered_procs;

    const int my_num_bodies = my_body_nsets.size();
    if (p == comm.MyPID())
    {
      for (int b = 0; b < my_num_bodies; ++b)
      {
        auto check = done.insert(b);
        if (check.second == false) continue;

        build_bid = b;

        std::vector<std::pair<int, int>> targets;

        connections.getOthers(b, targets, NULL);

        for (auto& t : targets)
        {
          sendto[t.second] = 1;
          already_considered_procs.insert(t.second);
        }

        mpi_requests.resize(targets.size());
        for (unsigned i = 0; i < targets.size(); ++i)
        {
          const int from_proc = p;
          const int to_proc = targets[i].second;
          const int* your_bid = &targets[i].first;
          MPI_Request& mpir = mpi_requests[i];

          exporter.ISend(from_proc, to_proc, your_bid, 1, 100, mpir);
        }

        // consider only one body at once
        break;
      }

      // increase the proc count as soon as all bodies of the current proc
      // have been considered
      if (build_bid == -1) ++next_proc;
    }

    // check if we can go to the next processor
    comm.Broadcast(&next_proc, 1, p);

    // go to the next processor
    if (p != next_proc)
    {
      p = next_proc;
      continue;
    }

    comm.SumAll(sendto.data(), receivefrom.data(), num_procs);
    int rsum = 0;
    for (int r : receivefrom) rsum += r;

    while (rsum > 0)
    {
      for (int i = 0; i < receivefrom[comm.MyPID()]; ++i)
      {
        std::vector<int> rdata;
        int rfrom_proc = -1;
        int rtag = -1;
        int rlength = -1;
        exporter.ReceiveAny(rfrom_proc, rtag, rdata, rlength);

        already_considered_procs.insert(rfrom_proc);

        if (rtag != 100) dserror("Wrong tag!");

        if (rlength != 1) dserror("Wrong length!");

        if (build_bid != -1 and build_bid != rdata[0]) dserror("Something went wrong!");

        build_bid = rdata[0];

        //        std::cout << "Proc #" << comm.MyPID() << " received " << rlength <<
        //            " node gid(s) from proc #" << rfrom_proc << std::endl;
      }

      if (build_bid != -1) done.insert(build_bid);

      // wait till all data were successfully sent and received
      for (auto& mpir : mpi_requests) exporter.Wait(mpir);

      std::vector<std::pair<int, int>> targets;
      connections.getOthers(build_bid, targets, &already_considered_procs);

      std::fill(sendto.data(), sendto.data() + num_procs, 0);
      for (auto& t : targets)
      {
        sendto[t.second] = 1;
        already_considered_procs.insert(t.second);
      }

      //      mpi_requests.clear();
      mpi_requests.resize(targets.size());
      for (unsigned i = 0; i < targets.size(); ++i)
      {
        const int from_proc = comm.MyPID();
        const int to_proc = targets[i].second;
        const int* your_bid = &targets[i].first;
        MPI_Request& mpir = mpi_requests[i];

        exporter.ISend(from_proc, to_proc, your_bid, 1, 100, mpir);
      }

      comm.SumAll(sendto.data(), receivefrom.data(), num_procs);

      // compute check sum
      rsum = 0;
      for (int r : receivefrom) rsum += r;
    }

    Teuchos::RCP<const node_set> ptr = Teuchos::null;
    if (build_bid >= 0)
      ptr = Teuchos::rcpFromRef(my_body_nsets[build_bid]);
    else
      ptr = Teuchos::rcp(new node_set());

    body_node_maps.push_back(createNodeMapFromBody(comm, *ptr, format));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::Connectivity::createNodeMapFromBody(
    const Epetra_Comm& comm, const node_set& my_node_set, const MapFormat format)
{
  const int mypid = comm.MyPID();

  std::vector<int> my_node_vec;
  my_node_vec.reserve(my_node_set.size());

  switch (format)
  {
    // the stored node sets are already in column format, just transform
    // the pointer set into an id set
    case MapFormat::column:
      std::transform(my_node_set.begin(), my_node_set.end(), std::back_inserter(my_node_vec),
          [](const DRT::Node* node) -> int { return node->Id(); });
      break;
    // ignore all ghosted nodes in this case
    case MapFormat::row:
      for (const DRT::Node* node : my_node_set)
        if (node->Owner() == mypid) my_node_vec.push_back(node->Id());
      break;
    default:
      dserror("Unknown MapFormat!");
      exit(EXIT_FAILURE);
  }

  // note: the ids are not yet correctly sorted, since the sets were based on
  // the pointer addresses rather than on the GIDs of the nodes
  std::sort(my_node_vec.begin(), my_node_vec.end());

  return Teuchos::rcp(new Epetra_Map(-1, my_node_vec.size(), my_node_vec.data(), 0, comm));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::Connectivity::receiveAnyOfMyGids(received_gid_set& all_received_ngids,
    const std::vector<node_set*>& my_new_body_nsets, const Epetra_Comm& comm)
{
  const int num_bodies = all_received_ngids.size();
  for (int i = 0; i < num_bodies; ++i)
  {
    auto& received_ngids = all_received_ngids[i];

    Teuchos::RCP<node_set> tmp_ptr = Teuchos::null;
    // point to the new node set
    if (my_new_body_nsets[i]) tmp_ptr = Teuchos::rcpFromRef(*my_new_body_nsets[i]);
    // create a new dummy set for NULL pointers
    else
      tmp_ptr = Teuchos::rcp(new node_set());
    auto& connected_nset = *tmp_ptr;

    const int mypid = comm.MyPID();
    const int num_procs = comm.NumProc();

    // find ghosted nodes
    std::unordered_map<int, std::vector<int>> my_ghosted_ngids;
    for (const DRT::Node* cnode : connected_nset)
    {
      const bool isghosted = cnode->Owner() != comm.MyPID();
      if (isghosted)
      {
        std::vector<int>& proc_ngids = my_ghosted_ngids[cnode->Owner()];
        proc_ngids.push_back(cnode->Id());
      }
    }

    // prepare the send operation:
    // communicate the number of sending procs to an individual proc
    std::vector<int> lrec(num_procs, 0);
    std::vector<int> grec(num_procs, 0);
    for (auto& pg : my_ghosted_ngids) ++lrec[pg.first];
    comm.SumAll(lrec.data(), grec.data(), num_procs);

    // build exporter
    DRT::Exporter exporter(comm);

    // send data
    std::vector<MPI_Request> mpi_requests(my_ghosted_ngids.size());
    auto mpir = mpi_requests.begin();
    for (auto pg = my_ghosted_ngids.begin(); pg != my_ghosted_ngids.end(); ++pg, ++mpir)
    {
      int to_proc = pg->first;
      int from_proc = mypid;

      exporter.ISend(from_proc, to_proc, pg->second.data(), pg->second.size(), 100, *mpir);
    }

    // receive data
    for (int i = 0; i < grec[mypid]; ++i)
    {
      std::vector<int> rdata;
      int rfrom_proc = -1;
      int rtag = -1;
      int rlength = -1;
      exporter.ReceiveAny(rfrom_proc, rtag, rdata, rlength);

      //      std::cout << "Proc #" << mypid << " received " << rlength << " node gids "
      //          "from proc #" << rfrom_proc << std::endl;

      if (rtag != 100) dserror("The received tag seems to be wrong! (tag=%d)", rtag);

      std::copy(rdata.begin(), rdata.end(), std::back_inserter(received_ngids[rfrom_proc]));
    }

    // wait till all data were successfully sent and received
    for (auto mpi_request : mpi_requests) exporter.Wait(mpi_request);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::Connectivity::fillMyConnectedNodeSet(
    node_set& connected_nset, const DRT::DiscretizationInterface& full_discret)
{
  std::set<int> done_element_ids;
  std::set<int> connected_node_gids;
  for (const DRT::Node* node : connected_nset) connected_node_gids.insert(node->Id());
  std::vector<const DRT::Node*> connected_nodes(connected_nset.begin(), connected_nset.end());

  unsigned i = 0;
  for (;;)
  {
    const DRT::Node* cnode = connected_nodes[i++];

    const DRT::Element* const* adj_eles = cnode->Elements();
    int num_adj_eles = cnode->NumElement();

    for (int e = 0; e < num_adj_eles; ++e)
    {
      const DRT::Element* ele = adj_eles[e];
      const auto echeck = done_element_ids.insert(ele->Id());
      if (echeck.second == false) continue;

      const DRT::Node* const* nodes = ele->Nodes();
      for (int n = 0; n < ele->NumNode(); ++n)
      {
        const DRT::Node* ele_node = nodes[n];
        const auto ncheck = connected_node_gids.insert(ele_node->Id());
        if (ncheck.second) connected_nodes.push_back(ele_node);
      }
    }

    // no new nodes have been added, thus we can stop here
    if (i == connected_nodes.size()) break;
  }

  connected_nset.insert(connected_nodes.begin(), connected_nodes.end());
}
