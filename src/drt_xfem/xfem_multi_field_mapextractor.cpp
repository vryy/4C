/*----------------------------------------------------------------------------*/
/** \file

\brief MultiFieldMapExtractor class to handle different discretizations
       with joint interfaces

\maintainer Matthias Mayr

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "xfem_multi_field_mapextractor.H"

#include "xfield_field_coupling_dofset.H"
#include "xfield_field_coupling.H"

#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_mapextractor.H"

#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_exporter.H"

#include "../drt_fsi/fsi_matrixtransform.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XFEM::MultiFieldMapExtractor::MultiFieldMapExtractor()
    : isinit_(false),
      issetup_(false),
      max_num_reserved_dofs_per_node_(0),
      comm_(Teuchos::null),
      slave_discret_vec_(0),
      element_map_extractor_(Teuchos::null),
      idiscret_(Teuchos::null),
      icoupl_dofset_(Teuchos::null)
{
  // intentionally left blank
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::Reset(unsigned num_dis, bool full)
{
  // set the flag to false (just to be sure)
  issetup_ = false;

  // --------------------------------------------------------------------------
  // reset the slave sided map extractors
  slave_map_extractors_.clear();
  slave_map_extractors_.resize(
      num_dis, std::vector<Teuchos::RCP<LINALG::MultiMapExtractor>>(NUM_MAP_TYPES, Teuchos::null));
  // loop over the number of discretizations
  for (unsigned i = 0; i < slave_map_extractors_.size(); ++i)
    // loop over the map types (0: DoF's, 1: Nodes)
    for (unsigned j = 0; j < NUM_MAP_TYPES; ++j)
      slave_map_extractors_[i][j] = Teuchos::rcp(new LINALG::MultiMapExtractor());

  // --------------------------------------------------------------------------
  // reset the master sided map extractor
  master_map_extractor_.clear();
  master_map_extractor_.resize(NUM_MAP_TYPES, Teuchos::null);
  // loop over the two map extractor types (0: DoF's, 1: Nodes)
  std::vector<Teuchos::RCP<LINALG::MultiMapExtractor>>::iterator it;
  for (it = master_map_extractor_.begin(); it != master_map_extractor_.end(); ++it)
    (*it) = Teuchos::RCP<LINALG::MultiMapExtractor>(new LINALG::MultiMapExtractor());

  // --------------------------------------------------------------------------
  // reset the interface coupling objects
  interface_couplings_.clear();
  interface_couplings_.resize(num_dis, Teuchos::null);
  std::vector<Teuchos::RCP<XFEM::XFieldField::Coupling>>::iterator iit;
  for (iit = interface_couplings_.begin(); iit != interface_couplings_.end(); ++iit)
    (*iit) = Teuchos::rcp(new XFEM::XFieldField::Coupling());

  // --------------------------------------------------------------------------
  // reset the element map extractor
  element_map_extractor_ = Teuchos::null;
  element_map_extractor_ = Teuchos::rcp<LINALG::MultiMapExtractor>(new LINALG::MultiMapExtractor());

  // clear these variables only if a full reset is desired!
  if (full)
  {
    // set the flag to false (just to be sure)
    isinit_ = false;

    xfem_dis_ids_.clear();

    slave_discret_vec_.clear();
    slave_discret_id_map_.clear();

    idiscret_ = Teuchos::null;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::Init(const XDisVec& dis_vec, int max_num_reserved_dofs_per_node)
{
  // reset flags
  isinit_ = false;
  issetup_ = false;

  // sanity check
  if (dis_vec.size() < 2)
    dserror(
        "You gave less than 2 discretizations. What are you planning "
        "to do? Seems as you are wrong here...");

  // get the communicator (supposed to be the same for all discretizations)
  comm_ = Teuchos::rcpFromRef(dis_vec[0]->Comm());

  max_num_reserved_dofs_per_node_ = max_num_reserved_dofs_per_node;

  // reset member variables
  Reset(dis_vec.size());

  // save the slave discretization vector
  slave_discret_vec_ = dis_vec;

  // build the slave discret id map
  BuildSlaveDiscretIdMap();

  // --------------------------------------------------------------------------
  // looking for xFEM discretizations
  // --------------------------------------------------------------------------
  XDisVec::const_iterator cit_dis;
  int d = 0;
  for (cit_dis = SlDisVec().begin(); cit_dis != SlDisVec().end(); ++cit_dis)
  {
    if (Teuchos::rcp_dynamic_cast<const DRT::DiscretizationXFEM>(*cit_dis).is_null())
      xfem_dis_ids_.insert(d++);
  }

  // --------------------------------------------------------------------------
  // get a std::set holding all interface node GIDs
  // (this has to be done only once, since it is independent of any
  //  redistribution)
  // --------------------------------------------------------------------------
  BuildGlobalInterfaceNodeGidSet();

  // ------------------------------------------------------------------------
  // create an auxiliary master interface discretization
  // ------------------------------------------------------------------------
  idiscret_ = Teuchos::rcp(
      new DRT::Discretization("multifield_interface", Teuchos::rcp<Epetra_Comm>(Comm().Clone())));

  // ------------------------------------------------------------------------
  // (1) create a list of coupling discretizations per node on this proc and
  //     pack it
  // ------------------------------------------------------------------------
  std::map<int, std::set<int>> my_coupled_sl_dis;
  std::map<int, std::set<int>> g_coupled_sl_dis;
  std::map<int, std::set<int>>::const_iterator cit_map;
  std::set<int>::const_iterator cit_set;

  for (std::set<int>::const_iterator ngid = GInterfaceNodeGidSet().begin();
       ngid != GInterfaceNodeGidSet().end(); ++ngid)
  {
    for (unsigned d = 0; d < SlDisVec().size(); ++d)
    {
      if (SlDisVec()[d]->NodeRowMap()->MyGID(*ngid))
      {
        std::pair<std::set<int>::iterator, bool> is_unique = my_coupled_sl_dis[*ngid].insert(d);
        if (not is_unique.second)
          dserror(
              "That is impossible, since this row node cannot be part of the "
              "same discretization on one single processor more than once!");
      }
    }
  }
  // collect the information over all processors
  std::vector<int> sendsize(2);
  std::vector<int> sendgid;
  std::vector<char> sendset;
  sendgid.reserve(my_coupled_sl_dis.size());

  std::vector<int> receivedsize(2);
  std::vector<int> receivedgid;
  std::vector<char> receivedset;

  DRT::PackBuffer data;
  // --------------------------------------------------------------------------
  // pack gids and std::set's separately
  // --------------------------------------------------------------------------
  // --- count set size
  for (cit_map = my_coupled_sl_dis.begin(); cit_map != my_coupled_sl_dis.end(); ++cit_map)
  {
    DRT::ParObject::AddtoPack(data, cit_map->second);
  }

  // --- activate packing
  data.StartPacking();

  // --- pack
  for (cit_map = my_coupled_sl_dis.begin(); cit_map != my_coupled_sl_dis.end(); ++cit_map)
  {
    sendgid.push_back(cit_map->first);
    DRT::ParObject::AddtoPack(data, cit_map->second);
  }

  // swap into std::vector<char>
  swap(sendset, data());

  // get the send information (how many messages will be received)
  sendsize[0] = sendgid.size();
  sendsize[1] = sendset.size();

  // ------------------------------------------------------------------------
  // (2) round robin loop
  // ------------------------------------------------------------------------
  // create an exporter for point to point communication
  DRT::Exporter exporter(Comm());
  const int numprocs = Comm().NumProc();

  for (int p = 0; p < numprocs; ++p)
  {
    // Send block to next proc. Receive a block from the last proc
    if (p > 0)
    {
      int myrank = Comm().MyPID();
      int tag = myrank;

      int frompid = myrank;
      int topid = (myrank + 1) % numprocs;

      MPI_Request sizerequest;
      exporter.ISend(frompid, topid, &sendsize[0], 2, tag, sizerequest);

      // send gid information
      MPI_Request gidrequest;
      exporter.ISend(frompid, topid, &sendgid[0], sendgid.size(), tag * 10, gidrequest);

      // send set data
      MPI_Request setrequest;
      exporter.ISend(frompid, topid, &sendset[0], sendset.size(), tag * 100, setrequest);

      // make sure that you do not think you received something if
      // you didn't
      if (not receivedset.empty() or not receivedgid.empty() or not receivedsize.empty())
        dserror("Received data objects are not empty!");

      // receive from predecessor
      frompid = (myrank + numprocs - 1) % numprocs;

      // receive size information
      int length = 0;
      exporter.ReceiveAny(frompid, tag, receivedsize, length);
      if (length != 2 or tag != frompid)
        dserror(
            "Size information got mixed up!\n"
            "Received length = %d, Expected length = %d \n"
            "Received tag    = %d, Expected tag    = %d",
            length, 2, tag, frompid);

      exporter.Wait(sizerequest);

      // receive the gids
      exporter.ReceiveAny(frompid, tag, receivedgid, length);
      if (length != receivedsize[0] or tag != frompid * 10)
        dserror(
            "GID information got mixed up! \n"
            "Received length = %d, Expected length = %d \n"
            "Received tag    = %d, Expected tag    = %d",
            length, receivedsize[0], tag, frompid * 10);

      exporter.Wait(gidrequest);

      // receive the sets
      exporter.ReceiveAny(frompid, tag, receivedset, length);
      if (length != receivedsize[1] or tag != frompid * 100)
        dserror(
            "Set information got mixed up! \n"
            "Received length = %d, Expected length = %d \n"
            "Received tag    = %d, Expected tag    = %d",
            length, receivedsize[1], tag, frompid * 100);

      exporter.Wait(setrequest);
    }
    // in the first step, we keep all nodes on the owning proc
    else
    {
      // no need to communicate in the first step
      swap(receivedsize, sendsize);
      swap(receivedgid, sendgid);
      swap(receivedset, sendset);
    }

    // ------------------------------------------------------------------------
    // (3) unpack received block
    // ------------------------------------------------------------------------
    std::vector<char>::size_type index = 0;
    int j = 0;
    while (index < receivedset.size())
    {
      int gid = receivedgid[j];
      // the set gets cleared at the beginning of the ExtractfromPack routine!
      std::set<int> rs;
      DRT::ParObject::ExtractfromPack(index, receivedset, rs);
      g_coupled_sl_dis[gid].insert(rs.begin(), rs.end());
      ++j;
    }

    // the received data will be sent to the next proc
    swap(receivedsize, sendsize);
    swap(receivedgid, sendgid);
    swap(receivedset, sendset);

    // we need a new receive buffer
    receivedsize.clear();
    receivedgid.clear();
    receivedset.clear();
  }

#if 0
  for (int p=0; p<numprocs; ++p)
  {
    if (Comm().MyPID()==p)
    {
      std::cout << "--------------- Processor " << p << " ---------------------" << std::endl;
      std::cout << "---------------  size = " << g_coupled_sl_dis.size() <<
          " ----------------------" << std::endl;
      for (cit_map=g_coupled_sl_dis.begin();cit_map!=g_coupled_sl_dis.end();++cit_map)
      {
        std::cout << "GID: " << cit_map->first << " | ";
        for (cit_set=cit_map->second.begin();cit_set!=cit_map->second.end();++cit_set)
          std::cout << *cit_set << "  ";
        std::cout << "\n";
      }
      std::cout << "\n\n\n\n";
    }
    Comm().Barrier();
  }
#endif

  // ID's for the master/slave coupling maps
  std::vector<std::vector<int>> my_master_inode_gids(SlDisVec().size(), std::vector<int>(0));

  // loop over all (global) coupling nodes
  for (cit_map = g_coupled_sl_dis.begin(); cit_map != g_coupled_sl_dis.end(); ++cit_map)
  {
    // current considered nodal GID
    int ngid = cit_map->first;

    // ------------------------------------------------------------------------
    // find the slave discretization ID we want to copy from
    // ------------------------------------------------------------------------
    int sl_dis_id_to_copy_from = -1;
    for (cit_set = cit_map->second.begin(); cit_set != cit_map->second.end(); ++cit_set)
    {
      // prefer and choose the first XFEM discretization (if there are more than
      // one joint at this node)
      if (IsXFemDis(*cit_set))
      {
        sl_dis_id_to_copy_from = *cit_set;
        break;
      }
    }
    /* If there is no XFEM discretization, we choose the first one. Since the
     * the coupled discretization set is ordered (by default), this should be
     * still deterministic. */
    if (sl_dis_id_to_copy_from == -1) sl_dis_id_to_copy_from = *cit_map->second.begin();

    // add the node on the owning processor
    if (SlDisVec()[sl_dis_id_to_copy_from]->NodeRowMap()->MyGID(ngid))
    {
      // clone the node, thus it becomes independent of any redistribution
      DRT::Node* node = SlDisVec()[sl_dis_id_to_copy_from]->gNode(ngid);
      Teuchos::RCP<DRT::Node> inode = Teuchos::rcp(node->Clone());
      idiscret_->AddNode(inode);
      // store the id for the master/slave coupling maps
      for (cit_set = cit_map->second.begin(); cit_set != cit_map->second.end(); ++cit_set)
        my_master_inode_gids[*cit_set].push_back(inode->Id());
    }
  }

  // build the master interface coupling node row maps
  BuildMasterInterfaceNodeMaps(my_master_inode_gids);

  // --------------------------------------------------------------------------
  // set the interface matrix transformation objects
  // --------------------------------------------------------------------------
  interface_matrix_row_transformers_.resize(NumSlDis(), Teuchos::null);
  interface_matrix_col_transformers_.resize(NumSlDis(), Teuchos::null);
  interface_matrix_row_col_transformers_.resize(NumSlDis(), Teuchos::null);
  for (unsigned i = 0; i < NumSlDis(); ++i)
  {
    interface_matrix_row_transformers_[i] = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform());
    interface_matrix_col_transformers_[i] = Teuchos::rcp(new FSI::UTILS::MatrixColTransform());
    interface_matrix_row_col_transformers_[i] =
        Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform());
  }
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::Setup()
{
  CheckInit();

  // reset the flag
  issetup_ = false;

  // first call reset (no full reset!)
  Reset(NumSlDis(), false);

  // build the slave node map extractor objects
  BuildSlaveNodeMapExtractors();

  // build the slave dof map extractor objects
  BuildSlaveDofMapExtractors();

  /* build the interface coupling dof set and set it in the auxiliary interface
   * discretization */
  BuildInterfaceCouplingDofSet();

  // build the master (i.e. auxiliary interface) node map extractor object
  BuildMasterNodeMapExtractor();

  // build the master (i.e. auxiliary interface) dof map extractor object
  BuildMasterDofMapExtractor();

  // build the interface coupling adapters
  BuildInterfaceCouplingAdapters();

  // build element map extractor
  BuildElementMapExtractor();

  // Everything is done. Set the flag.
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildSlaveNodeMapExtractors()
{
  XDisVec::const_iterator cit_dis;
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(2, Teuchos::null);

  // ------------------------------------------------------------------------
  /* build the interface row node GID maps of each slave discretization
   * on each single proc */
  // ------------------------------------------------------------------------
  unsigned dis_count = 0;
  for (cit_dis = SlDisVec().begin(); cit_dis != SlDisVec().end(); ++cit_dis)
  {
    std::vector<int> my_interface_row_node_gids(0);
    std::vector<int> my_non_interface_row_node_gids(0);

    const int num_my_rnodes = (*cit_dis)->NumMyRowNodes();
    int* my_row_node_gids = (*cit_dis)->NodeRowMap()->MyGlobalElements();

    for (unsigned nlid = 0; nlid < static_cast<unsigned>(num_my_rnodes); ++nlid)
    {
      int ngid = my_row_node_gids[nlid];

      // find the interface gids
      if (IsInterfaceNode(ngid))
        my_interface_row_node_gids.push_back(ngid);
      else
        my_non_interface_row_node_gids.push_back(ngid);
    }

    // slave sided interface node maps
    partial_maps[MULTIFIELD::block_interface] = Teuchos::null;
    partial_maps[MULTIFIELD::block_interface] =
        Teuchos::rcp(new Epetra_Map(-1, static_cast<int>(my_interface_row_node_gids.size()),
            &my_interface_row_node_gids[0], 0, Comm()));

    // slave sided non-interface node maps
    partial_maps[MULTIFIELD::block_non_interface] = Teuchos::null;
    partial_maps[MULTIFIELD::block_non_interface] =
        Teuchos::rcp(new Epetra_Map(-1, static_cast<int>(my_non_interface_row_node_gids.size()),
            &my_non_interface_row_node_gids[0], 0, Comm()));

    // setup node map extractor
    slave_map_extractors_[dis_count++][map_nodes]->Setup(*((*cit_dis)->NodeRowMap()), partial_maps);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildSlaveDofMapExtractors()
{
  std::vector<int> my_sl_interface_dofs(0);
  std::vector<int> my_sl_non_interface_dofs(0);

  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(2, Teuchos::null);

  // loop over all slave discretizations
  XDisVec::const_iterator cit_dis;
  unsigned dis_count = 0;
  for (cit_dis = SlDisVec().begin(); cit_dis != SlDisVec().end(); ++cit_dis)
  {
    int* my_node_gids = (*cit_dis)->NodeRowMap()->MyGlobalElements();

    // loop over my nodes
    for (int nlid = 0; nlid < (*cit_dis)->NodeRowMap()->NumMyElements(); ++nlid)
    {
      int ngid = my_node_gids[nlid];
      // ----------------------------------------------------------------------
      // interface DoF's
      // ----------------------------------------------------------------------
      if (SlaveNodeRowMap(dis_count, MULTIFIELD::block_interface).MyGID(ngid))
      {
        const DRT::Node* node = (*cit_dis)->lRowNode(nlid);
        const unsigned numdof = (*cit_dis)->NumDof(node);

        for (unsigned i = 0; i < numdof; ++i)
          my_sl_interface_dofs.push_back((*cit_dis)->Dof(node, i));
      }
      // ----------------------------------------------------------------------
      // non-interface DoF's
      // ----------------------------------------------------------------------
      else
      {
        const DRT::Node* node = (*cit_dis)->lRowNode(nlid);
        const unsigned numdof = (*cit_dis)->NumDof(node);

        for (unsigned i = 0; i < numdof; ++i)
          my_sl_non_interface_dofs.push_back((*cit_dis)->Dof(node, i));
      }
    }
    // create slave interface dof row map
    partial_maps[MULTIFIELD::block_interface] = Teuchos::null;
    partial_maps[MULTIFIELD::block_interface] = Teuchos::rcp(new Epetra_Map(
        -1, static_cast<int>(my_sl_interface_dofs.size()), &my_sl_interface_dofs[0], 0, Comm()));

    // create slave non-interface dof row map
    partial_maps[MULTIFIELD::block_non_interface] = Teuchos::null;
    partial_maps[MULTIFIELD::block_non_interface] =
        Teuchos::rcp(new Epetra_Map(-1, static_cast<int>(my_sl_non_interface_dofs.size()),
            &my_sl_non_interface_dofs[0], 0, Comm()));

    // setup dof map extractor
    slave_map_extractors_[dis_count++][map_dofs]->Setup(*((*cit_dis)->DofRowMap()), partial_maps);

    // clear the vectors for a new round
    my_sl_interface_dofs.clear();
    my_sl_non_interface_dofs.clear();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildInterfaceCouplingDofSet()
{
  std::map<int, int> sl_max_num_dof_per_inode;
  std::map<int, int> ma_max_num_dof_per_inode;
  int g_num_std_dof = -1;
  for (unsigned i = 0; i < NumSlDis(); ++i)
  {
    const Epetra_Map& sl_inodemap = SlaveNodeRowMap(i, MULTIFIELD::block_interface);
    const Epetra_Map& ma_inodemap = MasterInterfaceNodeRowMap(i);

    int nnodes = sl_inodemap.NumMyElements();
    int* ngids = sl_inodemap.MyGlobalElements();

    int my_num_std_dof = -1;

    for (int j = 0; j < nnodes; ++j)
    {
      const DRT::Node* node = SlDiscret(i).gNode(ngids[j]);
      const int numdof = SlDiscret(i).NumDof(node);
      my_num_std_dof = SlDiscret(i).NumStandardDof(0, node);
      sl_max_num_dof_per_inode[ngids[j]] = numdof;
    }

    DRT::Exporter export_max_dof_num(sl_inodemap, ma_inodemap, SlDiscret(i).Comm());
    export_max_dof_num.Export(sl_max_num_dof_per_inode);

    // communicate the number of standard DoF's
    // Supposed to be the same value on all discretizations and all nodes.
    if (g_num_std_dof == -1) Comm().MaxAll(&my_num_std_dof, &g_num_std_dof, 1);

    if (my_num_std_dof != -1 and g_num_std_dof != my_num_std_dof)
      dserror(
          "The number of standard DoF's is not equal on all procs"
          "and discretizations!");

    std::map<int, int>::const_iterator cit;
    for (cit = sl_max_num_dof_per_inode.begin(); cit != sl_max_num_dof_per_inode.end(); ++cit)
    {
      std::map<int, int>::iterator ma_pos = ma_max_num_dof_per_inode.find(cit->first);
      // If the nodal GID was not yet added, we insert the number of DoF's.
      if (ma_pos == ma_max_num_dof_per_inode.end())
      {
        ma_max_num_dof_per_inode[cit->first] = cit->second;
      }
      // If the nodal GID was already added, we calculate the maximum value.
      else
      {
        ma_pos->second = std::max(ma_pos->second, cit->second);
      }
    }

    // clear the slave maximum number of DoF's map for a new round
    sl_max_num_dof_per_inode.clear();
  }

  // get the node index range over all procs
  if (g_interface_node_gid_set_.empty()) dserror("set of global nodes is empty?!");
  const int min_all_gid = *(g_interface_node_gid_set_.begin());
  const int max_all_gid = *(g_interface_node_gid_set_.rbegin());

  const int g_node_index_range = max_all_gid - min_all_gid + 1;

  // create a new xfield/field coupling DoF set
  icoupl_dofset_ =
      Teuchos::rcp(new XFEM::XFieldField::CouplingDofSet(max_num_reserved_dofs_per_node_,
          g_node_index_range, g_num_std_dof, ma_max_num_dof_per_inode));

  // set the new dof-set and finish the interface discretization
  idiscret_->ReplaceDofSet(0, icoupl_dofset_, true);

  idiscret_->FillComplete(true, false, false);
#if (0)
  idiscret_->Print(std::cout);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildMasterNodeMapExtractor()
{
  Teuchos::RCP<Epetra_Map> fullmap = Teuchos::null;
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(2 * NumSlDis(), Teuchos::null);

  // --------------------------------------------------------------------------
  // interface DoF's
  // --------------------------------------------------------------------------
  for (unsigned i = 0; i < NumSlDis(); ++i) partial_maps.at(i) = master_interface_node_maps_[i];

  // --------------------------------------------------------------------------
  // non-interface DoF's
  // --------------------------------------------------------------------------
  for (unsigned i = 0; i < NumSlDis(); ++i)
    partial_maps.at(NumSlDis() + i) =
        SlMapExtractor(i, map_nodes).Map(MULTIFIELD::block_non_interface);

  // --------------------------------------------------------------------------
  // create non-overlapping full dof map
  // --------------------------------------------------------------------------
  fullmap = Teuchos::null;
  // add interface nodes
  fullmap = Teuchos::rcp<Epetra_Map>(new Epetra_Map(*IDiscret().NodeRowMap()));

  // merge non-interface nodes into the full map
  for (unsigned i = NumSlDis(); i < partial_maps.size(); ++i)
    fullmap = LINALG::MergeMap(*fullmap, *partial_maps[i], false);

  // setup map extractor
  master_map_extractor_[map_nodes]->Setup(*fullmap, partial_maps);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildMasterDofMapExtractor()
{
  std::vector<int> my_ma_interface_dofs(0);

  /* the 1-st num_dis partial maps are the master interface dof maps,
   * the 2-nd num_dis partial maps are the master non-interface dof maps
   * (which correspond to the slave_non_interface_dof_row_maps) */
  Teuchos::RCP<Epetra_Map> fullmap = Teuchos::null;
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(2 * NumSlDis(), Teuchos::null);

  // --------------------------------------------------------------------------
  // interface DoF's
  // --------------------------------------------------------------------------
  for (unsigned i = 0; i < NumSlDis(); ++i)
  {
    my_ma_interface_dofs.clear();

    const int num_my_inodes = MasterInterfaceNodeRowMap(i).NumMyElements();
    int* inode_gids = MasterInterfaceNodeRowMap(i).MyGlobalElements();

    // get the dofs of the master interface coupling nodes
    for (int nlid = 0; nlid < num_my_inodes; ++nlid)
    {
      int ngid = inode_gids[nlid];
      const DRT::Node* inode = IDiscret().gNode(ngid);
      const unsigned numdof = IDiscret().NumDof(inode);
      for (unsigned j = 0; j < numdof; ++j)
        my_ma_interface_dofs.push_back(IDiscret().Dof(inode, j));
    }
    partial_maps.at(i) = Teuchos::rcp<const Epetra_Map>(new Epetra_Map(
        -1, static_cast<int>(my_ma_interface_dofs.size()), &my_ma_interface_dofs[0], 0, Comm()));
  }

  // --------------------------------------------------------------------------
  // non-interface DoF's
  // --------------------------------------------------------------------------
  for (unsigned i = 0; i < NumSlDis(); ++i)
    partial_maps.at(NumSlDis() + i) =
        SlMapExtractor(i, map_dofs).Map(MULTIFIELD::block_non_interface);

  // --------------------------------------------------------------------------
  // create non-overlapping full dof map
  // --------------------------------------------------------------------------
  fullmap = Teuchos::null;

  // add interface DoF's
  fullmap = Teuchos::rcp<Epetra_Map>(new Epetra_Map(*IDiscret().DofRowMap()));

  // merge non-interface DoF's into the full map
  for (unsigned i = NumSlDis(); i < partial_maps.size(); ++i)
  {
    fullmap = LINALG::MergeMap(*fullmap, *partial_maps[i], false);
  }

  // setup map extractor
  master_map_extractor_[map_dofs]->Setup(*fullmap, partial_maps);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildInterfaceCouplingAdapters()
{
  CheckInit();

  std::vector<Teuchos::RCP<XFEM::XFieldField::Coupling>>::iterator it;
  unsigned i = 0;
  for (it = interface_couplings_.begin(); it != interface_couplings_.end(); ++it)
  {
    /* Set the slave discretization to the discretization with minimum number
     * of DoF's at each interface node. This is true by construction. */
    (*it)->Init(XFEM::XFieldField::Coupling::min_dof_slave);

    /* Setup the interface coupling objects. The interface discretization
     * is always the master discretization. Since the GID's at the interface
     * coincide in the coupling interface maps, the master interface map
     * becomes the permuted slave interface map. */
    (*it)->SetupCoupling(IDiscret(), SlDiscret(i), MasterInterfaceNodeRowMap(i),
        SlaveNodeRowMap(i, MULTIFIELD::block_interface), MasterInterfaceNodeRowMap(i), -1);

    ++i;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildElementMapExtractor()
{
  CheckInit();

  Teuchos::RCP<Epetra_Map> fullmap = Teuchos::null;
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(NumSlDis(), Teuchos::null);
  unsigned d = 0;
  XDisVec::const_iterator cit;
  for (cit = SlDisVec().begin(); cit != SlDisVec().end(); ++cit)
  {
    // get the element row map of each wrapped discretization
    partial_maps[d] = Teuchos::rcp((*cit)->ElementRowMap(), false);

    // merge the partial maps to the full map
    fullmap = LINALG::MergeMap(fullmap, partial_maps[d], false);

    // increase discretization counter
    ++d;
  }
  // setup the element multi map extractor
  element_map_extractor_->Setup(*fullmap, partial_maps);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::MultiFieldMapExtractor::ExtractVector(
    const Epetra_Vector& full, enum FieldName field, enum MapType map_type) const
{
  const int dis_id = SlaveId(field);

  /* the partial map is equivalent to the full map (of desired type) from
   * the field slave map extractor */
  const Teuchos::RCP<const Epetra_Map>& sl_full_map = SlMapExtractor(dis_id, map_type).FullMap();

  if (sl_full_map.is_null()) dserror("null full map for field %s", FieldName2String(field).c_str());

  // create a new vector
  Teuchos::RCP<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*sl_full_map));

  // extract the actual vector and return it
  ExtractVector(full, dis_id, *vec, map_type);

  return vec;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XFEM::MultiFieldMapExtractor::ExtractVector(
    const Epetra_MultiVector& full, enum FieldName field, enum MapType map_type) const
{
  const int dis_id = SlaveId(field);

  /* the partial map is equivalent to the full map (of desired type) from
   * the field slave map extractor */
  const Teuchos::RCP<const Epetra_Map>& sl_full_map = SlMapExtractor(dis_id, map_type).FullMap();

  if (sl_full_map.is_null()) dserror("null full map for field %s", FieldName2String(field).c_str());

  // create a new multi vector
  Teuchos::RCP<Epetra_MultiVector> vec =
      Teuchos::rcp(new Epetra_MultiVector(*sl_full_map, full.NumVectors()));

  // extract the actual vector and return it
  ExtractVector(full, dis_id, *vec, map_type);

  return vec;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::ExtractVector(const Epetra_MultiVector& full, int block,
    Epetra_MultiVector& partial, enum MapType map_type) const
{
  // --------------------------------------------------------------------------
  // extract the non-interface part
  // --------------------------------------------------------------------------
  Teuchos::RCP<Epetra_MultiVector> partial_non_interface =
      MaMapExtractor(map_type).ExtractVector(full, NumSlDis() + block);
  SlMapExtractor(block, map_type)
      .InsertVector(*partial_non_interface, MULTIFIELD::block_non_interface, partial);

  // --------------------------------------------------------------------------
  // extract the interface part
  // --------------------------------------------------------------------------
  Teuchos::RCP<Epetra_MultiVector> partial_ma_interface =
      MaMapExtractor(map_type).ExtractVector(full, block);
  Teuchos::RCP<Epetra_MultiVector> partial_sl_interface =
      ICoupling(block).MasterToSlave(partial_ma_interface, map_type);
  SlMapExtractor(block, map_type)
      .InsertVector(*partial_sl_interface, MULTIFIELD::block_interface, partial);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::ExtractElementVector(
    const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial) const
{
  element_map_extractor_->ExtractVector(full, block, partial);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::InsertElementVector(
    const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full) const
{
  element_map_extractor_->InsertVector(partial, block, full);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::AddElementVector(
    const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full, double scale) const
{
  element_map_extractor_->AddVector(partial, block, full, scale);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::MultiFieldMapExtractor::InsertVector(
    const Epetra_Vector& partial, enum FieldName field, enum MapType map_type) const
{
  const int dis_id = SlaveId(field);
  Teuchos::RCP<Epetra_Vector> full = Teuchos::rcp(new Epetra_Vector(*FullMap(map_type)));
  InsertVector(partial, dis_id, *full, map_type);
  return full;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XFEM::MultiFieldMapExtractor::InsertVector(
    const Epetra_MultiVector& partial, enum FieldName field, enum MapType map_type) const
{
  const int dis_id = SlaveId(field);

  Teuchos::RCP<Epetra_MultiVector> full =
      Teuchos::rcp(new Epetra_MultiVector(*FullMap(map_type), partial.NumVectors()));

  InsertVector(partial, dis_id, *full, map_type);
  return full;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::InsertVector(const Epetra_MultiVector& partial, int block,
    Epetra_MultiVector& full, enum MapType map_type) const
{
  // --------------------------------------------------------------------------
  // insert the non_interface part
  // --------------------------------------------------------------------------
  Teuchos::RCP<Epetra_MultiVector> partial_non_interface =
      SlMapExtractor(block, map_type).ExtractVector(partial, MULTIFIELD::block_non_interface);

  MaMapExtractor(map_type).InsertVector(*partial_non_interface, NumSlDis() + block, full);

  // --------------------------------------------------------------------------
  // insert the interface part
  // --------------------------------------------------------------------------
  Teuchos::RCP<Epetra_MultiVector> partial_sl_interface =
      SlMapExtractor(block, map_type).ExtractVector(partial, MULTIFIELD::block_interface);

  Teuchos::RCP<Epetra_MultiVector> partial_ma_interface =
      ICoupling(block).SlaveToMaster(partial_sl_interface, map_type);

  MaMapExtractor(map_type).InsertVector(*partial_ma_interface, block, full);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::AddVector(const Epetra_MultiVector& partial, int block,
    Epetra_MultiVector& full, double scale, enum MapType map_type) const
{
  // --------------------------------------------------------------------------
  // insert the non_interface part
  // --------------------------------------------------------------------------
  Teuchos::RCP<Epetra_MultiVector> partial_non_interface =
      SlMapExtractor(block, map_type).ExtractVector(partial, MULTIFIELD::block_non_interface);

  MaMapExtractor(map_type).AddVector(*partial_non_interface, NumSlDis() + block, full, scale);

  // --------------------------------------------------------------------------
  // insert the interface part
  // --------------------------------------------------------------------------
  Teuchos::RCP<Epetra_MultiVector> partial_sl_interface =
      SlMapExtractor(block, map_type).ExtractVector(partial, MULTIFIELD::block_interface);

  Teuchos::RCP<Epetra_MultiVector> partial_ma_interface =
      ICoupling(block).SlaveToMaster(partial_sl_interface, map_type);

  MaMapExtractor(map_type).AddVector(*partial_ma_interface, block, full, scale);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildSlaveDiscretIdMap()
{
  XDisVec::const_iterator cit_dis;
  int dis_count = 0;
  for (cit_dis = SlDisVec().begin(); cit_dis != SlDisVec().end(); ++cit_dis)
  {
    const std::string& name = (*cit_dis)->Name();
    if (boost::iequals(name, "structure"))
      slave_discret_id_map_[structure] = dis_count;
    else if (boost::iequals(name, "xstructure"))
      slave_discret_id_map_[xstructure] = dis_count;
    else
      dserror("Unknown field discretization name \"%s\"!", name.c_str());
    // increase counter
    ++dis_count;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XFEM::MultiFieldMapExtractor::IsInterfaceNode(const int& ngid) const
{
  return (GInterfaceNodeGidSet().find(ngid) != GInterfaceNodeGidSet().end());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XFEM::MultiFieldMapExtractor::IsXFemDis(int dis_id) const
{
  return (xfem_dis_ids_.find(dis_id) != xfem_dis_ids_.end());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::MultiFieldMapExtractor::SlaveId(enum FieldName field) const
{
  std::map<enum FieldName, int>::const_iterator cit = slave_discret_id_map_.find(field);
  if (cit == slave_discret_id_map_.end())
    dserror("The slave field \"%s\" could not be found!", FieldName2String(field).c_str());

  return cit->second;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildGlobalInterfaceNodeGidSet()
{
  Comm().Barrier();
  std::set<int> g_unique_row_node_gid_set;
  g_interface_node_gid_set_.clear();

  // loop over all proc's
  for (unsigned p = 0; p < static_cast<unsigned>(Comm().NumProc()); ++p)
  {
    int num_my_unique_row_nodes = 0;
    std::vector<int> my_unique_row_node_gid_vec(0);

    int num_my_interface_row_nodes = 0;
    std::vector<int> my_interface_row_node_gid_vec(0);

    if (p == static_cast<unsigned>(Comm().MyPID()))
    {
      std::set<int> my_unique_row_node_gid_set;
      std::set<int> my_interface_row_node_gid_set;
      for (unsigned j = 0; j < slave_discret_vec_.size(); ++j)
      {
        for (unsigned i = 0; i < static_cast<unsigned>(slave_discret_vec_[j]->NumMyRowNodes()); ++i)
        {
          int gid = slave_discret_vec_[j]->NodeRowMap()->GID(i);
          // insert the gid and check if it is unique on the current processor
          std::pair<std::set<int>::iterator, bool> is_unique =
              my_unique_row_node_gid_set.insert(gid);
          if (not is_unique.second) my_interface_row_node_gid_set.insert(gid);
        }
      }
      // copy the unique set into a vector
      num_my_unique_row_nodes = my_unique_row_node_gid_set.size();
      my_unique_row_node_gid_vec.resize(num_my_unique_row_nodes);

      std::copy(my_unique_row_node_gid_set.begin(), my_unique_row_node_gid_set.end(),
          my_unique_row_node_gid_vec.begin());

      // copy the interface set into a vector
      num_my_interface_row_nodes = my_interface_row_node_gid_set.size();
      my_interface_row_node_gid_vec.resize(num_my_interface_row_nodes);

      std::copy(my_interface_row_node_gid_set.begin(), my_interface_row_node_gid_set.end(),
          my_interface_row_node_gid_vec.begin());
    }
    // wait since only one proc did all the work
    Comm().Barrier();
    // ------------------------------------------------------------------------
    // send the unique row node GID vector from processor p to all proc's
    // ------------------------------------------------------------------------
    Comm().Broadcast(&num_my_unique_row_nodes, 1, p);
    if (num_my_unique_row_nodes == 0) continue;
    my_unique_row_node_gid_vec.resize(num_my_unique_row_nodes, -1);
    Comm().Broadcast(&my_unique_row_node_gid_vec[0], num_my_unique_row_nodes, p);

    // ------------------------------------------------------------------------
    // send the interface row node GID vector from processor p to all proc's
    // ------------------------------------------------------------------------
    Comm().Broadcast(&num_my_interface_row_nodes, 1, p);
    if (num_my_interface_row_nodes > 0)
    {
      my_interface_row_node_gid_vec.resize(num_my_interface_row_nodes, -1);
      Comm().Broadcast(&my_interface_row_node_gid_vec[0], num_my_interface_row_nodes, p);
      // create/extend the global interface row node gid set
      g_interface_node_gid_set_.insert(
          my_interface_row_node_gid_vec.begin(), my_interface_row_node_gid_vec.end());
    }
    // ------------------------------------------------------------------------
    /* Insert the local unique row node GIDs into a global set. If the GID was
     * already inserted by another proc, we have to extend our global interface
     * row node set as well. */
    // ------------------------------------------------------------------------
    for (std::vector<int>::const_iterator cit = my_unique_row_node_gid_vec.begin();
         cit != my_unique_row_node_gid_vec.end(); ++cit)
    {
      std::pair<std::set<int>::iterator, bool> is_unique = g_unique_row_node_gid_set.insert(*cit);
      if (not is_unique.second) g_interface_node_gid_set_.insert(*cit);
    }
  }  // end: loop over all proc's
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::BuildMasterInterfaceNodeMaps(
    const std::vector<std::vector<int>>& my_master_interface_node_gids)
{
  for (unsigned i = 0; i < my_master_interface_node_gids.size(); ++i)
  {
    master_interface_node_maps_.push_back(
        Teuchos::rcp(new Epetra_Map(-1, static_cast<int>(my_master_interface_node_gids[i].size()),
            &my_master_interface_node_gids[i][0], 0, Comm())));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> XFEM::MultiFieldMapExtractor::NodeRowMap(
    enum FieldName field, enum MULTIFIELD::BlockType block) const
{
  switch (block)
  {
    case MULTIFIELD::block_interface:
    {
      return Teuchos::rcpFromRef<const Epetra_Map>(MasterInterfaceNodeRowMap(field));
      break;
    }
    case MULTIFIELD::block_non_interface:
    {
      return Teuchos::rcpFromRef<const Epetra_Map>(
          SlaveNodeRowMap(field, MULTIFIELD::block_non_interface));
      break;
    }
    default:
      dserror("Unknown block type!");
      exit(EXIT_FAILURE);
  }
  // hoops, shouldn't happen ...
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::AddMatrix(const LINALG::SparseOperator& partial_mat, int block,
    LINALG::SparseOperator& full_mat, double scale)
{
  const LINALG::BlockSparseMatrixBase* block_mat =
      dynamic_cast<const LINALG::BlockSparseMatrixBase*>(&partial_mat);
  if (not block_mat) dserror("The partial matrix must be a  LINALG::BlockSparseMatrix!");
  if (block_mat->Rows() != 2 or block_mat->Cols() != 2)
    dserror("We support only 2x2 block matrices!");

  LINALG::SparseMatrix* sp_mat = dynamic_cast<LINALG::SparseMatrix*>(&full_mat);
  if (not sp_mat) dserror("The full matrix must be a LINALG::SparseMatrix!");

  AddMatrix(*block_mat, block, *sp_mat, scale);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::AddMatrix(const LINALG::BlockSparseMatrixBase& partial_mat,
    int block, LINALG::SparseMatrix& full_mat, double scale)
{
  CheckInitSetup();
  // --------------------------------------------------------------------------
  // non-interface DoF's
  // --------------------------------------------------------------------------
  // Add block non_interface/non_interface. Here is no communication necessary.
  /* ToDo Maybe there is a way to circumvent the add command by using some kind
   *      of direct assignment for the non-interface DoF's? The current
   *      implementation makes it necessary to allocate almost the double
   *      amount of memory!                                        hiermeier */
  full_mat.Add(partial_mat.Matrix(MULTIFIELD::block_non_interface, MULTIFIELD::block_non_interface),
      false, scale, 1.0);

  // --------------------------------------------------------------------------
  // interface DoF's
  // --------------------------------------------------------------------------
  // (0) Add block non_interface/interface
  const LINALG::SparseMatrix& src_ni =
      partial_mat.Matrix(MULTIFIELD::block_non_interface, MULTIFIELD::block_interface);
  IMatColTransform(block)(partial_mat.FullRowMap(), partial_mat.FullColMap(), src_ni, scale,
      ADAPTER::CouplingSlaveConverter(ICoupling(block)), full_mat, false, true);

  // (1) Add block interface/non_interface
  const LINALG::SparseMatrix& src_in =
      partial_mat.Matrix(MULTIFIELD::block_interface, MULTIFIELD::block_non_interface);
  IMatRowTransform(block)(
      src_in, scale, ADAPTER::CouplingSlaveConverter(ICoupling(block)), full_mat, true);

  // (2) Add block interface/interface
  const LINALG::SparseMatrix& src_ii =
      partial_mat.Matrix(MULTIFIELD::block_interface, MULTIFIELD::block_interface);
  IMatRowColTransform(block)(src_ii, scale, ADAPTER::CouplingSlaveConverter(ICoupling(block)),
      ADAPTER::CouplingSlaveConverter(ICoupling(block)), full_mat, false, true);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& XFEM::MultiFieldMapExtractor::FullMap(
    enum MapType map_type) const
{
  return MaMapExtractor(map_type).FullMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Node* XFEM::MultiFieldMapExtractor::gINode(const int& gid) const
{
  return IDiscret().gNode(gid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map* XFEM::MultiFieldMapExtractor::INodeRowMap() const
{
  return IDiscret().NodeRowMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::MultiFieldMapExtractor::INumDof(const DRT::Node* inode) const
{
  return IDiscret().NumDof(0, inode);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::MultiFieldMapExtractor::INumStandardDof() const
{
  return icoupl_dofset_->NumStandardDofPerNode();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::MultiFieldMapExtractor::IDof(const DRT::Node* inode, int dof) const
{
  return IDiscret().Dof(0, inode, dof);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::IDof(const DRT::Node* inode, std::vector<int>& dofs) const
{
  IDiscret().Dof(static_cast<unsigned>(0), inode, dofs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::MultiFieldMapExtractor::IDof(std::vector<int>& dof, DRT::Node* inode,
    unsigned nodaldofset_id, const DRT::Element* element) const
{
  IDiscret().Dof(dof, inode, 0, nodaldofset_id, element);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& XFEM::MultiFieldMapExtractor::SlaveNodeRowMap(
    unsigned dis_id, enum MULTIFIELD::BlockType btype) const
{
  CheckInit();
  return *(SlMapExtractor(dis_id, map_nodes).Map(btype));
}
