/*---------------------------------------------------------------------*/
/*! \file
\brief Routines for extending contact interface ghosting by Round-Robin loop

\level 2


*/
/*---------------------------------------------------------------------*/

#include "contact_interface.H"

#include "contact_element.H"
#include "friction_node.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

/*----------------------------------------------------------------------*
 |  store the required ghosting within a round              farah 10/13 |
 |  robin iteration for current interface (public)                      |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::RoundRobinExtendGhosting(bool firstevaluation)
{
  std::vector<int> element_GIDs_to_be_ghosted;
  std::vector<int> node_GIDs_to_be_ghosted;
  for (int k = 0; k < SlaveColElements()->NumMyElements(); ++k)
  {
    int gid = SlaveColElements()->GID(k);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i", gid);
    CoElement* slave_ele = dynamic_cast<CoElement*>(ele);

    for (int j = 0; j < slave_ele->MoData().NumSearchElements(); ++j)
    {
      int gid2 = slave_ele->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %", gid2);
      CoElement* melement = dynamic_cast<CoElement*>(ele2);

      element_GIDs_to_be_ghosted.push_back(melement->Id());

      for (int z = 0; z < melement->NumNode(); ++z)
      {
        int gidn = melement->NodeIds()[z];
        node_GIDs_to_be_ghosted.push_back(gidn);
      }
    }
    // reset found elements
    slave_ele->DeleteSearchElements();
  }

  Teuchos::RCP<Epetra_Map> currently_ghosted_elements = Teuchos::rcp(new Epetra_Map(
      -1, (int)element_GIDs_to_be_ghosted.size(), &element_GIDs_to_be_ghosted[0], 0, Comm()));
  Teuchos::RCP<Epetra_Map> currently_ghosted_nodes = Teuchos::rcp(new Epetra_Map(
      -1, (int)node_GIDs_to_be_ghosted.size(), &node_GIDs_to_be_ghosted[0], 0, Comm()));

  if (firstevaluation)
  {
    eextendedghosting_ = Teuchos::rcp(new Epetra_Map(*currently_ghosted_elements));
    nextendedghosting_ = Teuchos::rcp(new Epetra_Map(*currently_ghosted_nodes));
  }
  else
  {
    eextendedghosting_ = LINALG::MergeMap(eextendedghosting_, currently_ghosted_elements, true);
    nextendedghosting_ = LINALG::MergeMap(nextendedghosting_, currently_ghosted_nodes, true);
  }

  return;
}

/*----------------------------------------------------------------------*
 | perform the ownership change within a round robin        farah 10/13 |
 | iteration                                                            |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::RoundRobinChangeOwnership()
{
  // pack/unpack friction nodes is only required for wear problems
  // so we should create a slightly redundant function for the wear-
  // interface and exclude the friction node packing from here

  // get friction type
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");

  // change master-side proc ownership
  // some local variables
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(Comm().Clone());
  const int myrank = comm->MyPID();
  const int numproc = comm->NumProc();
  const int torank = (myrank + 1) % numproc;              // to
  const int fromrank = (myrank + numproc - 1) % numproc;  // from

  // new elements, new nodes
  std::vector<int> ncol, nrow;
  std::vector<int> ecol, erow;

  // create dummy
  Teuchos::RCP<Epetra_Map> MasterColNodesdummy = Teuchos::rcp(new Epetra_Map(*MasterColNodes()));
  Teuchos::RCP<Epetra_Map> MasterColelesdummy = Teuchos::rcp(new Epetra_Map(*MasterColElements()));

  // create origin maps
  Teuchos::RCP<Epetra_Map> SCN = Teuchos::rcp(new Epetra_Map(*SlaveColNodes()));
  Teuchos::RCP<Epetra_Map> SCE = Teuchos::rcp(new Epetra_Map(*SlaveColElements()));
  Teuchos::RCP<Epetra_Map> SRN = Teuchos::rcp(new Epetra_Map(*SlaveRowNodes()));
  Teuchos::RCP<Epetra_Map> SRE = Teuchos::rcp(new Epetra_Map(*SlaveRowElements()));

  // *****************************************
  // Elements
  // *****************************************
  std::vector<char> sdataeles;
  std::vector<char> rdataeles;

  // vector containing all proc ids
  std::vector<int> allproc(numproc);
  for (int i = 0; i < numproc; ++i) allproc[i] = i;

  // get exporter
  DRT::Exporter exporter(idiscret_->Comm());

  // create data buffer
  DRT::PackBuffer dataeles;

  // pack data - first just reserving the mem.
  for (int i = 0; i < (int)MasterColelesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColelesdummy->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i", gid);
    MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

    mele->Pack(dataeles);

    int ghost = 0;
    DRT::ParObject::AddtoPack(dataeles, ghost);
  }

  dataeles.StartPacking();

  // now pack/store
  for (int i = 0; i < (int)MasterColelesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColelesdummy->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i", gid);
    MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

    mele->Pack(dataeles);

    // check for ghosting
    int ghost;
    if (mele->Owner() == myrank)
      ghost = 1;
    else
      ghost = 0;

    DRT::ParObject::AddtoPack(dataeles, ghost);
  }
  std::swap(sdataeles, dataeles());

  // delete the elements from discretization
  for (int i = 0; i < (int)MasterColelesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColelesdummy->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i", gid);
    MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

    // check for ghosting
    if (mele->Owner() == myrank)
    {
      idiscret_->DeleteElement(mele->Id());
    }
  }

  // send the information
  MPI_Request request;
  exporter.ISend(myrank, torank, &(sdataeles[0]), (int)sdataeles.size(), 1234, request);

  // receive the information
  int length = rdataeles.size();
  int tag = -1;
  int from = -1;
  exporter.ReceiveAny(from, tag, rdataeles, length);
  if (tag != 1234 or from != fromrank)
    dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank,
        from, myrank);

  // ---- unpack ----
  {
    // Put received nodes into discretization
    std::vector<char>::size_type index = 0;
    while (index < rdataeles.size())
    {
      std::vector<char> data;
      int ghost = -1;
      DRT::ParObject::ExtractfromPack(index, rdataeles, data);
      DRT::ParObject::ExtractfromPack(index, rdataeles, ghost);
      if (ghost == -1) dserror("UNPACK ERROR!!!!!!!!!");

      // this Teuchos::rcp holds the memory of the ele
      Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data), true);
      Teuchos::RCP<MORTAR::MortarElement> ele =
          Teuchos::rcp_dynamic_cast<MORTAR::MortarElement>(object);
      if (ele == Teuchos::null) dserror("Received object is not an ele");

      // add whether its a row ele
      if (ghost == 1)
      {
        ele->SetOwner(myrank);
        idiscret_->AddElement(ele);

        // to new ele
        erow.push_back(ele->Id());
        ecol.push_back(ele->Id());
      }
      else
      {
        ecol.push_back(ele->Id());
      }
    }
  }  // end unpack

  // wait for all communication to finish
  exporter.Wait(request);
  comm->Barrier();

  // *****************************************
  // NODES
  // *****************************************
  std::vector<char> sdatanodes;
  std::vector<char> rdatanodes;

  // get exporter
  DRT::Exporter exportern(idiscret_->Comm());

  DRT::PackBuffer datanodes;

  // pack data -- col map --> should prevent further ghosting!
  for (int i = 0; i < (int)MasterColNodesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColNodesdummy->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find ele with gid %i", gid);

    if (ftype == INPAR::CONTACT::friction_none)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(node);
      cnode->Pack(datanodes);
    }
    else
    {
      FriNode* cnode = dynamic_cast<FriNode*>(node);
      cnode->Pack(datanodes);
    }

    int ghost = 0;
    DRT::ParObject::AddtoPack(datanodes, ghost);
  }

  datanodes.StartPacking();
  for (int i = 0; i < (int)MasterColNodesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColNodesdummy->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find ele with gid %i", gid);

    // check for ghosting
    int ghost;

    if (ftype == INPAR::CONTACT::friction_none)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(node);
      cnode->Pack(datanodes);

      if (cnode->Owner() == myrank)
        ghost = 1;
      else
        ghost = 0;
    }
    else
    {
      FriNode* cnode = dynamic_cast<FriNode*>(node);
      cnode->Pack(datanodes);

      if (cnode->Owner() == myrank)
        ghost = 1;
      else
        ghost = 0;
    }

    DRT::ParObject::AddtoPack(datanodes, ghost);
  }
  std::swap(sdatanodes, datanodes());

  // DELETING
  for (int i = 0; i < (int)MasterColNodesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColNodesdummy->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find ele with gid %i", gid);

    if (ftype == INPAR::CONTACT::friction_none)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(node);
      if (cnode->Owner() == myrank) idiscret_->DeleteNode(cnode->Id());
    }
    else
    {
      FriNode* cnode = dynamic_cast<FriNode*>(node);
      if (cnode->Owner() == myrank) idiscret_->DeleteNode(cnode->Id());
    }
  }

  // ---- send ----
  MPI_Request requestn;
  exportern.ISend(myrank, torank, &(sdatanodes[0]), (int)sdatanodes.size(), 1234, requestn);

  // ---- receive ----
  int lengthn = rdatanodes.size();
  int tagn = -1;
  int fromn = -1;
  exportern.ReceiveAny(fromn, tagn, rdatanodes, lengthn);
  if (tagn != 1234 or fromn != fromrank)
    dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank,
        fromn, myrank);

  // ---- unpack ----
  {
    // Put received nodes into discretization
    std::vector<char>::size_type index = 0;
    while (index < rdatanodes.size())
    {
      std::vector<char> data;

      int ghost = -1;
      DRT::ParObject::ExtractfromPack(index, rdatanodes, data);
      DRT::ParObject::ExtractfromPack(index, rdatanodes, ghost);
      if (ghost == -1) dserror("UNPACK ERROR!!!!!!!!!");

      // this Teuchos::rcp holds the memory of the node
      Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data), true);

      if (ftype == INPAR::CONTACT::friction_none)
      {
        Teuchos::RCP<MORTAR::MortarNode> node =
            Teuchos::rcp_dynamic_cast<MORTAR::MortarNode>(object);
        if (node == Teuchos::null) dserror("Received object is not a node");

        if (ghost == 1)
        {
          node->SetOwner(myrank);
          idiscret_->AddNode(node);

          nrow.push_back(node->Id());
          ncol.push_back(node->Id());
        }
        else
        {
          // all others to col
          ncol.push_back(node->Id());
        }
      }
      else  // if friction...
      {
        Teuchos::RCP<FriNode> node = Teuchos::rcp_dynamic_cast<FriNode>(object);
        if (node == Teuchos::null) dserror("Received object is not a node");

        if (ghost == 1)
        {
          node->SetOwner(myrank);
          idiscret_->AddNode(node);

          nrow.push_back(node->Id());
          ncol.push_back(node->Id());
        }
        else
        {
          // all others to col
          ncol.push_back(node->Id());
        }
      }
    }
  }  // end unpack

  // wait for all communication to finish
  exportern.Wait(requestn);
  comm->Barrier();

  // create maps from sending
  Teuchos::RCP<Epetra_Map> noderowmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)nrow.size(), &nrow[0], 0, Comm()));
  Teuchos::RCP<Epetra_Map> nodecolmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)ncol.size(), &ncol[0], 0, Comm()));

  Teuchos::RCP<Epetra_Map> elerowmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)erow.size(), &erow[0], 0, Comm()));
  Teuchos::RCP<Epetra_Map> elecolmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)ecol.size(), &ecol[0], 0, Comm()));

  // Merge s/m column maps for eles and nodes
  Teuchos::RCP<Epetra_Map> colnodesfull = LINALG::MergeMap(nodecolmap, SCN, true);
  Teuchos::RCP<Epetra_Map> colelesfull = LINALG::MergeMap(elecolmap, SCE, true);

  // Merge s/m row maps for eles and nodes
  Teuchos::RCP<Epetra_Map> rownodesfull = LINALG::MergeMap(noderowmap, SRN, false);
  Teuchos::RCP<Epetra_Map> rowelesfull = LINALG::MergeMap(elerowmap, SRE, false);

  // to discretization
  // export nodes and elements to the row map
  Discret().ExportRowNodes(*rownodesfull);
  Discret().ExportRowElements(*rowelesfull);

  // export nodes and elements to the col map
  Discret().ExportColumnNodes(*colnodesfull);
  Discret().ExportColumnElements(*colelesfull);

  // ********************************************
  // call the (very) expensive FILLCOMPLETE()!
  // ********************************************
  // make sure discretization is complete
  FillComplete(true);

  return;
}


/*----------------------------------------------------------------------*
 |  change master ownership clockwise for contact            farah 10/13|
 |  interface without evaluation of the interface                       |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::RoundRobinDetectGhosting()
{
  if (SearchAlg() == INPAR::MORTAR::search_bfele)
    EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg() == INPAR::MORTAR::search_binarytree)
    EvaluateSearchBinarytree();
  else
    dserror("ERROR: Invalid search algorithm");

  // first ghosting for std. distribution
  RoundRobinExtendGhosting(true);

  // Init Maps
  Teuchos::RCP<Epetra_Map> initial_slave_node_column_map =
      Teuchos::rcp(new Epetra_Map(*SlaveColNodes()));
  Teuchos::RCP<Epetra_Map> initial_slave_element_column_map =
      Teuchos::rcp(new Epetra_Map(*SlaveColElements()));
  Teuchos::RCP<Epetra_Map> initial_master_node_column_map =
      Teuchos::rcp(new Epetra_Map(*MasterColNodes()));
  Teuchos::RCP<Epetra_Map> initial_master_element_column_map =
      Teuchos::rcp(new Epetra_Map(*MasterColElements()));

  // *************************************
  // start RR loop for current interface
  // *************************************
  // loop over all procs
  if (Comm().NumProc() > 1)
    for (int proc = 0; proc < (int)(Comm().NumProc()); ++proc)
    {
      // status output
      if (Comm().MyPID() == 0 && proc == 0) std::cout << "Round-Robin-Iteration #" << proc;
      if (Comm().MyPID() == 0 && proc > 0) std::cout << " #" << proc;

      // perform the ownership change
      RoundRobinChangeOwnership();

      // build new search tree or do nothing for bruteforce
      if (SearchAlg() == INPAR::MORTAR::search_binarytree)
        CreateSearchTree();
      else if (SearchAlg() != INPAR::MORTAR::search_bfele)
        dserror("ERROR: Invalid search algorithm");

      // evaluate interfaces
      if (proc < (int)(Comm().NumProc() - 1))
      {
        if (SearchAlg() == INPAR::MORTAR::search_bfele)
          EvaluateSearchBruteForce(SearchParam());
        else if (SearchAlg() == INPAR::MORTAR::search_binarytree)
          EvaluateSearchBinarytree();
        else
          dserror("ERROR: Invalid search algorithm");

        // other ghostings per iteration
        RoundRobinExtendGhosting(false);
      }
      else
      {
        // do nothing -- just switch
      }
    }  // end RR

  // Append ghost nodes/elements to be ghosted to existing column maps
  eextendedghosting_ = LINALG::MergeMap(eextendedghosting_, initial_slave_element_column_map, true);
  nextendedghosting_ = LINALG::MergeMap(nextendedghosting_, initial_slave_node_column_map, true);
  eextendedghosting_ =
      LINALG::MergeMap(eextendedghosting_, initial_master_element_column_map, true);
  nextendedghosting_ = LINALG::MergeMap(nextendedghosting_, initial_master_node_column_map, true);

  // finally extend ghosting
  Discret().ExportColumnElements(*eextendedghosting_);
  Discret().ExportColumnNodes(*nextendedghosting_);
  FillComplete(true);

  // reset extended ghosting maps
  eextendedghosting_ = Teuchos::null;
  nextendedghosting_ = Teuchos::null;

  // build new search tree or do nothing for bruteforce
  if (SearchAlg() == INPAR::MORTAR::search_binarytree)
    CreateSearchTree();
  else if (SearchAlg() != INPAR::MORTAR::search_bfele)
    dserror("ERROR: Invalid search algorithm");

  // final output for loop
  if (Comm().MyPID() == 0) std::cout << " Round-Robin loop done!" << std::endl;

  return;
}
