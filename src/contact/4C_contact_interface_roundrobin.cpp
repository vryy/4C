/*---------------------------------------------------------------------*/
/*! \file
\brief Routines for extending contact interface ghosting by Round-Robin loop

\level 2


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_interface.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  store the required ghosting within a round              farah 10/13 |
 |  robin iteration for current interface (public)                      |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::round_robin_extend_ghosting(bool firstevaluation)
{
  std::vector<int> element_GIDs_to_be_ghosted;
  std::vector<int> node_GIDs_to_be_ghosted;
  for (int k = 0; k < slave_col_elements()->NumMyElements(); ++k)
  {
    int gid = slave_col_elements()->GID(k);
    Core::Elements::Element* ele = discret().g_element(gid);
    if (!ele) FOUR_C_THROW("Cannot find ele with gid %i", gid);
    Element* slave_ele = dynamic_cast<Element*>(ele);

    for (int j = 0; j < slave_ele->mo_data().num_search_elements(); ++j)
    {
      int gid2 = slave_ele->mo_data().search_elements()[j];
      Core::Elements::Element* ele2 = idiscret_->g_element(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      element_GIDs_to_be_ghosted.push_back(melement->id());

      for (int z = 0; z < melement->num_node(); ++z)
      {
        int gidn = melement->node_ids()[z];
        node_GIDs_to_be_ghosted.push_back(gidn);
      }
    }
    // reset found elements
    slave_ele->delete_search_elements();
  }

  Teuchos::RCP<Epetra_Map> currently_ghosted_elements = Teuchos::rcp(new Epetra_Map(-1,
      (int)element_GIDs_to_be_ghosted.size(), element_GIDs_to_be_ghosted.data(), 0, get_comm()));
  Teuchos::RCP<Epetra_Map> currently_ghosted_nodes = Teuchos::rcp(new Epetra_Map(
      -1, (int)node_GIDs_to_be_ghosted.size(), node_GIDs_to_be_ghosted.data(), 0, get_comm()));

  if (firstevaluation)
  {
    eextendedghosting_ = Teuchos::rcp(new Epetra_Map(*currently_ghosted_elements));
    nextendedghosting_ = Teuchos::rcp(new Epetra_Map(*currently_ghosted_nodes));
  }
  else
  {
    eextendedghosting_ =
        Core::LinAlg::MergeMap(eextendedghosting_, currently_ghosted_elements, true);
    nextendedghosting_ = Core::LinAlg::MergeMap(nextendedghosting_, currently_ghosted_nodes, true);
  }

  return;
}

/*----------------------------------------------------------------------*
 | perform the ownership change within a round robin        farah 10/13 |
 | iteration                                                            |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::round_robin_change_ownership()
{
  // pack/unpack friction nodes is only required for wear problems
  // so we should create a slightly redundant function for the wear-
  // interface and exclude the friction node packing from here

  // get friction type
  Inpar::CONTACT::FrictionType ftype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(interface_params(), "FRICTION");

  // change master-side proc ownership
  // some local variables
  Teuchos::RCP<Epetra_Comm> comm_v = Teuchos::rcp(get_comm().Clone());
  const int myrank = comm_v->MyPID();
  const int numproc = comm_v->NumProc();
  const int torank = (myrank + 1) % numproc;              // to
  const int fromrank = (myrank + numproc - 1) % numproc;  // from

  // new elements, new nodes
  std::vector<int> ncol, nrow;
  std::vector<int> ecol, erow;

  // create dummy
  Teuchos::RCP<Epetra_Map> MasterColNodesdummy = Teuchos::rcp(new Epetra_Map(*master_col_nodes()));
  Teuchos::RCP<Epetra_Map> MasterColelesdummy =
      Teuchos::rcp(new Epetra_Map(*master_col_elements()));

  // create origin maps
  Teuchos::RCP<Epetra_Map> SCN = Teuchos::rcp(new Epetra_Map(*slave_col_nodes()));
  Teuchos::RCP<Epetra_Map> SCE = Teuchos::rcp(new Epetra_Map(*slave_col_elements()));
  Teuchos::RCP<Epetra_Map> SRN = Teuchos::rcp(new Epetra_Map(*slave_row_nodes()));
  Teuchos::RCP<Epetra_Map> SRE = Teuchos::rcp(new Epetra_Map(*slave_row_elements()));

  // *****************************************
  // Elements
  // *****************************************
  std::vector<char> sdataeles;
  std::vector<char> rdataeles;

  // vector containing all proc ids
  std::vector<int> allproc(numproc);
  for (int i = 0; i < numproc; ++i) allproc[i] = i;

  // get exporter
  Core::Communication::Exporter exporter(idiscret_->get_comm());

  // create data buffer
  Core::Communication::PackBuffer dataeles;

  // now pack/store
  for (int i = 0; i < (int)MasterColelesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColelesdummy->GID(i);
    Core::Elements::Element* ele = discret().g_element(gid);
    if (!ele) FOUR_C_THROW("Cannot find ele with gid %i", gid);
    Mortar::Element* mele = dynamic_cast<Mortar::Element*>(ele);

    mele->pack(dataeles);

    // check for ghosting
    const int ghost = (mele->owner() == myrank) ? 1 : 0;
    Core::Communication::ParObject::add_to_pack(dataeles, ghost);
  }
  std::swap(sdataeles, dataeles());

  // delete the elements from discretization
  for (int i = 0; i < (int)MasterColelesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColelesdummy->GID(i);
    Core::Elements::Element* ele = discret().g_element(gid);
    if (!ele) FOUR_C_THROW("Cannot find ele with gid %i", gid);
    Mortar::Element* mele = dynamic_cast<Mortar::Element*>(ele);

    // check for ghosting
    if (mele->owner() == myrank)
    {
      idiscret_->delete_element(mele->id());
    }
  }

  // send the information
  MPI_Request request;
  exporter.i_send(myrank, torank, sdataeles.data(), (int)sdataeles.size(), 1234, request);

  // receive the information
  int length = rdataeles.size();
  int tag = -1;
  int from = -1;
  exporter.receive_any(from, tag, rdataeles, length);
  if (tag != 1234 or from != fromrank)
    FOUR_C_THROW("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank,
        from, myrank);

  // ---- unpack ----
  {
    // Put received nodes into discretization
    std::vector<char>::size_type index = 0;
    while (index < rdataeles.size())
    {
      std::vector<char> data;
      int ghost = -1;
      Core::Communication::ParObject::extract_from_pack(index, rdataeles, data);
      Core::Communication::ParObject::extract_from_pack(index, rdataeles, ghost);
      if (ghost == -1) FOUR_C_THROW("UNPACK ERROR!!!!!!!!!");

      // this Teuchos::rcp holds the memory of the ele
      Teuchos::RCP<Core::Communication::ParObject> object =
          Teuchos::rcp(Core::Communication::Factory(data), true);
      Teuchos::RCP<Mortar::Element> ele = Teuchos::rcp_dynamic_cast<Mortar::Element>(object);
      if (ele == Teuchos::null) FOUR_C_THROW("Received object is not an ele");

      // add whether its a row ele
      if (ghost == 1)
      {
        ele->set_owner(myrank);
        idiscret_->add_element(ele);

        // to new ele
        erow.push_back(ele->id());
        ecol.push_back(ele->id());
      }
      else
      {
        ecol.push_back(ele->id());
      }
    }
  }  // end unpack

  // wait for all communication to finish
  exporter.wait(request);
  comm_v->Barrier();

  // *****************************************
  // NODES
  // *****************************************
  std::vector<char> sdatanodes;
  std::vector<char> rdatanodes;

  // get exporter
  Core::Communication::Exporter exportern(idiscret_->get_comm());

  Core::Communication::PackBuffer datanodes;

  for (int i = 0; i < (int)MasterColNodesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColNodesdummy->GID(i);
    Core::Nodes::Node* node = discret().g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find ele with gid %i", gid);

    // check for ghosting
    int ghost;

    if (ftype == Inpar::CONTACT::friction_none)
    {
      Mortar::Node* cnode = dynamic_cast<Mortar::Node*>(node);
      cnode->pack(datanodes);

      if (cnode->owner() == myrank)
        ghost = 1;
      else
        ghost = 0;
    }
    else
    {
      FriNode* cnode = dynamic_cast<FriNode*>(node);
      cnode->pack(datanodes);

      if (cnode->owner() == myrank)
        ghost = 1;
      else
        ghost = 0;
    }

    Core::Communication::ParObject::add_to_pack(datanodes, ghost);
  }
  std::swap(sdatanodes, datanodes());

  // DELETING
  for (int i = 0; i < (int)MasterColNodesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColNodesdummy->GID(i);
    Core::Nodes::Node* node = discret().g_node(gid);
    if (!node) FOUR_C_THROW("Cannot find ele with gid %i", gid);

    if (ftype == Inpar::CONTACT::friction_none)
    {
      Mortar::Node* cnode = dynamic_cast<Mortar::Node*>(node);
      if (cnode->owner() == myrank) idiscret_->delete_node(cnode->id());
    }
    else
    {
      FriNode* cnode = dynamic_cast<FriNode*>(node);
      if (cnode->owner() == myrank) idiscret_->delete_node(cnode->id());
    }
  }

  // ---- send ----
  MPI_Request requestn;
  exportern.i_send(myrank, torank, sdatanodes.data(), (int)sdatanodes.size(), 1234, requestn);

  // ---- receive ----
  int lengthn = rdatanodes.size();
  int tagn = -1;
  int fromn = -1;
  exportern.receive_any(fromn, tagn, rdatanodes, lengthn);
  if (tagn != 1234 or fromn != fromrank)
    FOUR_C_THROW("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank,
        fromn, myrank);

  // ---- unpack ----
  {
    // Put received nodes into discretization
    std::vector<char>::size_type index = 0;
    while (index < rdatanodes.size())
    {
      std::vector<char> data;

      int ghost = -1;
      Core::Communication::ParObject::extract_from_pack(index, rdatanodes, data);
      Core::Communication::ParObject::extract_from_pack(index, rdatanodes, ghost);
      if (ghost == -1) FOUR_C_THROW("UNPACK ERROR!!!!!!!!!");

      // this Teuchos::rcp holds the memory of the node
      Teuchos::RCP<Core::Communication::ParObject> object =
          Teuchos::rcp(Core::Communication::Factory(data), true);

      if (ftype == Inpar::CONTACT::friction_none)
      {
        Teuchos::RCP<Mortar::Node> node = Teuchos::rcp_dynamic_cast<Mortar::Node>(object);
        if (node == Teuchos::null) FOUR_C_THROW("Received object is not a node");

        if (ghost == 1)
        {
          node->set_owner(myrank);
          idiscret_->add_node(node);

          nrow.push_back(node->id());
          ncol.push_back(node->id());
        }
        else
        {
          // all others to col
          ncol.push_back(node->id());
        }
      }
      else  // if friction...
      {
        Teuchos::RCP<FriNode> node = Teuchos::rcp_dynamic_cast<FriNode>(object);
        if (node == Teuchos::null) FOUR_C_THROW("Received object is not a node");

        if (ghost == 1)
        {
          node->set_owner(myrank);
          idiscret_->add_node(node);

          nrow.push_back(node->id());
          ncol.push_back(node->id());
        }
        else
        {
          // all others to col
          ncol.push_back(node->id());
        }
      }
    }
  }  // end unpack

  // wait for all communication to finish
  exportern.wait(requestn);
  comm_v->Barrier();

  // create maps from sending
  Teuchos::RCP<Epetra_Map> noderowmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)nrow.size(), nrow.data(), 0, get_comm()));
  Teuchos::RCP<Epetra_Map> nodecolmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)ncol.size(), ncol.data(), 0, get_comm()));

  Teuchos::RCP<Epetra_Map> elerowmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)erow.size(), erow.data(), 0, get_comm()));
  Teuchos::RCP<Epetra_Map> elecolmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)ecol.size(), ecol.data(), 0, get_comm()));

  // Merge s/m column maps for eles and nodes
  Teuchos::RCP<Epetra_Map> colnodesfull = Core::LinAlg::MergeMap(nodecolmap, SCN, true);
  Teuchos::RCP<Epetra_Map> colelesfull = Core::LinAlg::MergeMap(elecolmap, SCE, true);

  // Merge s/m row maps for eles and nodes
  Teuchos::RCP<Epetra_Map> rownodesfull = Core::LinAlg::MergeMap(noderowmap, SRN, false);
  Teuchos::RCP<Epetra_Map> rowelesfull = Core::LinAlg::MergeMap(elerowmap, SRE, false);

  // to discretization
  // export nodes and elements to the row map
  discret().export_row_nodes(*rownodesfull);
  discret().export_row_elements(*rowelesfull);

  // export nodes and elements to the col map
  discret().export_column_nodes(*colnodesfull);
  discret().export_column_elements(*colelesfull);

  // ********************************************
  // call the (very) expensive FILLCOMPLETE()!
  // ********************************************
  // make sure discretization is complete
  fill_complete(true);

  return;
}


/*----------------------------------------------------------------------*
 |  change master ownership clockwise for contact            farah 10/13|
 |  interface without evaluation of the interface                       |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::round_robin_detect_ghosting()
{
  if (search_alg() == Inpar::Mortar::search_bfele)
    evaluate_search_brute_force(search_param());
  else if (search_alg() == Inpar::Mortar::search_binarytree)
    evaluate_search_binarytree();
  else
    FOUR_C_THROW("Invalid search algorithm");

  // first ghosting for std. distribution
  round_robin_extend_ghosting(true);

  // Init Maps
  Teuchos::RCP<Epetra_Map> initial_slave_node_column_map =
      Teuchos::rcp(new Epetra_Map(*slave_col_nodes()));
  Teuchos::RCP<Epetra_Map> initial_slave_element_column_map =
      Teuchos::rcp(new Epetra_Map(*slave_col_elements()));
  Teuchos::RCP<Epetra_Map> initial_master_node_column_map =
      Teuchos::rcp(new Epetra_Map(*master_col_nodes()));
  Teuchos::RCP<Epetra_Map> initial_master_element_column_map =
      Teuchos::rcp(new Epetra_Map(*master_col_elements()));

  // *************************************
  // start RR loop for current interface
  // *************************************
  // loop over all procs
  if (get_comm().NumProc() > 1)
    for (int proc = 0; proc < (int)(get_comm().NumProc()); ++proc)
    {
      // status output
      if (get_comm().MyPID() == 0 && proc == 0) std::cout << "Round-Robin-Iteration #" << proc;
      if (get_comm().MyPID() == 0 && proc > 0) std::cout << " #" << proc;

      // perform the ownership change
      round_robin_change_ownership();

      // build new search tree or do nothing for bruteforce
      if (search_alg() == Inpar::Mortar::search_binarytree)
        create_search_tree();
      else if (search_alg() != Inpar::Mortar::search_bfele)
        FOUR_C_THROW("Invalid search algorithm");

      // evaluate interfaces
      if (proc < (int)(get_comm().NumProc() - 1))
      {
        if (search_alg() == Inpar::Mortar::search_bfele)
          evaluate_search_brute_force(search_param());
        else if (search_alg() == Inpar::Mortar::search_binarytree)
          evaluate_search_binarytree();
        else
          FOUR_C_THROW("Invalid search algorithm");

        // other ghostings per iteration
        round_robin_extend_ghosting(false);
      }
      else
      {
        // do nothing -- just switch
      }
    }  // end RR

  // Append ghost nodes/elements to be ghosted to existing column maps
  eextendedghosting_ =
      Core::LinAlg::MergeMap(eextendedghosting_, initial_slave_element_column_map, true);
  nextendedghosting_ =
      Core::LinAlg::MergeMap(nextendedghosting_, initial_slave_node_column_map, true);
  eextendedghosting_ =
      Core::LinAlg::MergeMap(eextendedghosting_, initial_master_element_column_map, true);
  nextendedghosting_ =
      Core::LinAlg::MergeMap(nextendedghosting_, initial_master_node_column_map, true);

  // finally extend ghosting
  discret().export_column_elements(*eextendedghosting_);
  discret().export_column_nodes(*nextendedghosting_);
  fill_complete(true);

  // reset extended ghosting maps
  eextendedghosting_ = Teuchos::null;
  nextendedghosting_ = Teuchos::null;

  // build new search tree or do nothing for bruteforce
  if (search_alg() == Inpar::Mortar::search_binarytree)
    create_search_tree();
  else if (search_alg() != Inpar::Mortar::search_bfele)
    FOUR_C_THROW("Invalid search algorithm");

  // final output for loop
  if (get_comm().MyPID() == 0) std::cout << " Round-Robin loop done!" << std::endl;

  return;
}

FOUR_C_NAMESPACE_CLOSE
