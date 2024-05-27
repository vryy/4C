/*----------------------------------------------------------------------*/
/*! \file
\brief Search tree for unbiased self-contact problems

\level 2

*/
/*-----------------------------------------------------------------------*/

#include "4C_contact_selfcontact_binarytree_unbiased.hpp"

#include "4C_contact_element.hpp"
#include "4C_contact_node.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_utils_reference_configuration.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor UnbiasedSelfBinaryTree (public)                   schmidt 01/19|
 *----------------------------------------------------------------------*/
CONTACT::UnbiasedSelfBinaryTree::UnbiasedSelfBinaryTree(DRT::Discretization& discret,
    const Teuchos::ParameterList& iparams, Teuchos::RCP<Epetra_Map> elements, int dim, double eps)
    : SelfBinaryTree(discret, iparams, elements, dim, eps),
      two_half_pass_(iparams.get<bool>("Two_half_pass")),
      check_nonsmooth_selfcontactsurface_(iparams.get<bool>("Check_nonsmooth_selfcontactsurface")),
      searchele_all_proc_(iparams.get<bool>("Searchele_AllProc"))
{
  // safety check
  if (!two_half_pass_) FOUR_C_THROW("Only implemented for the two half pass approach so far!");

  return;
}

/*----------------------------------------------------------------------*
 | Add tree nodes to the contact pairs vector  (private)   schmidt 01/19|
 *----------------------------------------------------------------------*/
void CONTACT::UnbiasedSelfBinaryTree::add_tree_nodes_to_contact_pairs(
    Teuchos::RCP<SelfBinaryTreeNode> treenode1, Teuchos::RCP<SelfBinaryTreeNode> treenode2)
{
  bool addcontactpair(true);

  // check reference configuration for non smooth self contact
  if (two_half_pass_ and check_nonsmooth_selfcontactsurface_)
    addcontactpair = rough_check_ref_config(treenode1->elelist()[0], treenode2->elelist()[0]);

  if (addcontactpair)
  {
    set_contact_pairs()[treenode1->elelist()[0]].push_back(treenode2->elelist()[0]);
    set_contact_pairs()[treenode2->elelist()[0]].push_back(treenode1->elelist()[0]);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Calculate dual graph processor specific  (private)      schmidt 01/19|
 *----------------------------------------------------------------------*/
void CONTACT::UnbiasedSelfBinaryTree::calculate_proc_specific_dual_graph(
    std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>* dualGraph,
    const std::vector<int>& elelist, const int p)
{
  // loop over all self contact elements
  for (unsigned i = 0; i < elelist.size(); ++i)
  {
    // global id of current element
    const int gid = elelist[i];

    // initialize
    // ... vector of adjacent tree nodes (elements) of current element
    // ... vector of adjacent dual edges containing current element
    // ... vector of global IDs of possible adjacent elements
    std::vector<Teuchos::RCP<SelfBinaryTreeNode>> adjtreenodes;
    std::vector<Teuchos::RCP<SelfDualEdge>> adjdualedges;
    std::vector<int> possadjids;

    // get current elements and its nodes
    DRT::Element* element = discret().gElement(gid);
    if (!element) FOUR_C_THROW("Cannot find element with gid %\n", gid);
    DRT::Node** nodes = element->Nodes();
    if (!nodes) FOUR_C_THROW("Null pointer!");

    // skip further steps, if element is not owned by processor p
    if (element->Owner() != p) continue;

    // first tree node of one dual edge which includes current element
    // is the element itself saved as a tree node
    Teuchos::RCP<SelfBinaryTreeNode> node1 = leafsmap()[gid];

    // for 2D only: get the finite element nodes of the element and save them as end nodes of the
    // tree node
    if (dim() == 2)
    {
      std::vector<int> nodeIds;
      nodeIds.push_back(element->NodeIds()[0]);
      nodeIds.push_back(element->NodeIds()[1]);
      node1->SetEndnodes(nodeIds);
    }

    // get element specific first order nodes
    const int numnode = get_ele_specific_num_nodes(element);

    // loop over all first-order nodes of current element (here we make use of the fact that
    // first-order nodes are always stored before higher-order nodes)
    for (int j = 0; j < numnode; ++j)
    {
      DRT::Node* node = nodes[j];
      if (!node) FOUR_C_THROW("Null pointer!");

      // adjacent elements of current node
      int numE = node->NumElement();
      DRT::Element** adjElements = node->Elements();
      if (!adjElements) FOUR_C_THROW("Null pointer!");

      // loop over all adjacent elements of current node
      for (int k = 0; k < numE; ++k)
      {
        // get k-th adjacent element
        DRT::Element* adjElementk = adjElements[k];

        // we only need to collect information if current adjacent element is owned by processor p
        if (adjElementk->Owner() != p) continue;

        calculate_adjacent_tree_nodes_and_dual_edges(
            possadjids, gid, adjElementk, node1, adjtreenodes, adjdualedges);
      }  // all adjacent elements
    }    // all nodes

    // add the vector of adjacent tree nodes to the adjacency matrix we only need the matrix in 3D,
    // because in 2D the adjacency test works by comparing end nodes only
    if (dim() == 3) set_adjacencymatrix()[gid] = adjtreenodes;

    // get adjacent dual edges
    for (unsigned k = 0; k < adjdualedges.size(); ++k)
      for (unsigned j = 0; j < adjdualedges.size(); ++j)
        if (j != k) (*dualGraph)[adjdualedges[k]].push_back(adjdualedges[j]);
  }  // all elements

  return;
}  // calculate_proc_specific_dual_graph

/*----------------------------------------------------------------------*
 | Define search elements based on contact pairs (private) schmidt 01/19|
 *----------------------------------------------------------------------*/
void CONTACT::UnbiasedSelfBinaryTree::define_search_elements()
{
  const int eleID = contact_pairs().begin()->first;

  // do as long as there are still contact pairs
  if (contact_pairs().find(eleID) != contact_pairs().end() && !contact_pairs().empty())
  {
    // get the current element to content of "isslave"
    DRT::Element* element = discret().gElement(eleID);
    CONTACT::Element* celement = dynamic_cast<CONTACT::Element*>(element);
    if (celement->IsSlave() != true)
      FOUR_C_THROW("Element: this should not happen!");
    else
      for (int i = 0; i < element->num_node(); ++i)
      {
        DRT::Node* node = element->Nodes()[i];
        CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
        if (cnode->IsSlave() != true) FOUR_C_THROW("Node: this should not happen!");
      }

    // get the ID of elements in contact with current one
    std::vector<int> contacteleID = contact_pairs()[eleID];

    // erase the current element from list of contact pairs
    set_contact_pairs().erase(eleID);

    for (unsigned j = 0; j < contacteleID.size(); ++j)
    {
      celement->AddSearchElements(contacteleID[j]);

      // recursively call this function again
      if (contact_pairs().find(contacteleID[j]) != contact_pairs().end()) define_search_elements();
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Get the contracted node (private)                      schmidt 01/19|
 *----------------------------------------------------------------------*/
void CONTACT::UnbiasedSelfBinaryTree::get_contracted_node(
    Teuchos::RCP<SelfDualEdge>& contractedEdge, Teuchos::RCP<SelfBinaryTreeNode>& contractedNode)
{
  // call the base class
  CONTACT::SelfBinaryTree::get_contracted_node(contractedEdge, contractedNode);

  // add owner of contracted node
  contractedNode->SetParentOwner(
      contractedEdge->GetNode1()->Owner(), contractedEdge->GetNode2()->Owner());

  return;
}

/*----------------------------------------------------------------------*
 | initialize the unbiased self binary tree  (public)      schmidt 01/19|
 *----------------------------------------------------------------------*/
void CONTACT::UnbiasedSelfBinaryTree::Init()
{
  // call initialization method of base class
  MORTAR::BaseBinaryTree::Init();

  // initialize internal variables
  init_internal_variables();

  // calculate min. element length and set enlargement accordingly
  set_enlarge();

  // initialize binary tree leaf nodes and create element list
  std::vector<int> elelist;
  init_leaf_nodes_and_map(elelist);

  // initialize and calculate processor specific dual graph
  std::map<int, std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>>
      procdualgraph;

  // loop over all interface processors
  for (int p = 0; p < comm().NumProc(); ++p)
  {
    std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>> dualgraph;

    calculate_proc_specific_dual_graph(&dualgraph, elelist, p);
    procdualgraph[p] = dualgraph;
  }
  // now initialize unbiased self binary tree in a bottom-up way based on the processor specific
  // dual graph
  initialize_tree_bottom_up(&procdualgraph);

  return;
}

/*----------------------------------------------------------------------*
 | Initialize tree bottom-up based on processor specific                |
 | dual graph (private)                                    schmidt 01/19|
 *----------------------------------------------------------------------*/
void CONTACT::UnbiasedSelfBinaryTree::initialize_tree_bottom_up(
    std::map<int, std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>>*
        procdualGraph)
{
  // vector collecting root nodes
  set_roots().resize(0);

  std::map<int,
      std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>>::iterator
      iterator = procdualGraph->begin();
  std::map<int,
      std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>>::iterator
      iterator_end = procdualGraph->end();

  while (iterator != iterator_end)
  {
    std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>> dualGraph =
        iterator->second;

    //**********************************************************************
    // the idea is to empty the dual graph step by step
    //**********************************************************************
    while (!dualGraph.empty())
    {
      // get the edge with lowest costs (= the first edge in the dual graph as
      // the map is automatically sorted by the costs)  to contract it
      std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>::iterator iter =
          dualGraph.begin();
      Teuchos::RCP<SelfDualEdge> contractedEdge = iter->first;
      Teuchos::RCP<SelfBinaryTreeNode> newNode(Teuchos::null);
      get_contracted_node(contractedEdge, newNode);

      // update dualGraph
      // this means we have to create new edges, which include the new tree node and delete the
      // edges, which are adjacent to the contracted edge and update their neighbors additionally we
      // have to check if the tree is nearly complete

      // get the adjacent edges of the contracted edge
      std::vector<Teuchos::RCP<SelfDualEdge>> adjEdges = iter->second;

      // check if the new tree node includes whole self contact-surface in this case the tree node
      // has saved itself as adjacent edge
      if (adjEdges[0] == contractedEdge)
      {
        // save the tree node as root and continue the loop
        set_roots().push_back(newNode);
        dualGraph.erase(contractedEdge);
        continue;
      }

      update_dual_graph(contractedEdge, adjEdges, newNode, &dualGraph);
    }  // while(!(dualGraph).empty())
    //**********************************************************************
    ++iterator;
  }

  // complete the tree starting from its roots (top-down)
  if (roots().size() == 0) FOUR_C_THROW("No root tree node found!");

  // SAFETY CHECK: check whether all leaf nodes were assigned a root node
  std::vector<int> list = roots()[0]->elelist();
  for (unsigned i = 1; i < roots().size(); ++i)
  {
    std::vector<int> listi = roots()[i]->elelist();
    for (unsigned j = 0; j < listi.size(); ++j) list.push_back(listi[j]);
  }
  if (list.size() != leafsmap().size())
    FOUR_C_THROW(
        "Not all leaf nodes got assigned to a root node! This is not acceptable, but might appear "
        "for the unbiased self contact tree as only surface elements (nodes) are contracted (and "
        "thereby assigned to a root node in this algorithm), if they are owned by the same "
        "processor. So in the rare case when there are <= 2 surface elements that are connected "
        "and owned by one processor appear at one of the contact surfaces, those elements can not "
        "be assigned to a root node and the algorithm fails! %i elements were not assigned to a "
        "root node!\n\n"
        "POSSIBLE SOLUTION: choose a different number of processors!",
        leafsmap().size() - list.size());

  for (unsigned k = 0; k < roots().size(); ++k) roots()[k]->CompleteTree(0, enlarge());
  // output to screen
  if (comm().MyPID() == 0)
    std::cout << "\nFound " << roots().size() << " root node(s) for unbiased self binary tree."
              << std::endl;

  // in 3D we have to calculate adjacent tree nodes
  if (dim() == 3)
  {
    calculate_adjacent_leaves();
    calculate_adjacent_tnodes();
  }

  return;
}

/*----------------------------------------------------------------------*
 | Checks roughly whether self contact of two elements                  |
 | shall be evaluated  (private)                           schmidt 10/18|
 *----------------------------------------------------------------------*/
bool CONTACT::UnbiasedSelfBinaryTree::rough_check_ref_config(int ele1gid, int ele2gid)
{
  // variables
  bool refconfig(false);
  static CORE::LINALG::Matrix<2, 1> xicele1(true);
  static CORE::LINALG::Matrix<3, 1> ele1coords(true);
  static CORE::LINALG::Matrix<3, 1> ele1normal(true);
  static CORE::LINALG::Matrix<2, 1> xicele2(true);
  static CORE::LINALG::Matrix<3, 1> ele2coords(true);
  static CORE::LINALG::Matrix<3, 1> ele1ele2vec(true);
  static CORE::LINALG::Matrix<1, 1> scalarprod(true);

  // get center and normal of leaf1-element
  const DRT::Element* ele1 = discret().gElement(ele1gid);
  const CORE::FE::CellType dtele1 = ele1->Shape();
  switch (dtele1)
  {
    case CORE::FE::CellType::tri3:
    {
      xicele1.PutScalar(1.0 / 3.0);
      DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::tri3>(
          ele1, xicele1, ele1coords);
      DRT::UTILS::ComputeUnitNormalAtXiRefConfig<CORE::FE::CellType::tri3>(
          ele1, xicele1, ele1normal);
    }
    break;
    case CORE::FE::CellType::tri6:
    {
      xicele1.PutScalar(1.0 / 3.0);
      DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::tri6>(
          ele1, xicele1, ele1coords);
      DRT::UTILS::ComputeUnitNormalAtXiRefConfig<CORE::FE::CellType::tri6>(
          ele1, xicele1, ele1normal);
    }
    break;
    case CORE::FE::CellType::quad4:
    {
      xicele1.PutScalar(0.0);
      DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::quad4>(
          ele1, xicele1, ele1coords);
      DRT::UTILS::ComputeUnitNormalAtXiRefConfig<CORE::FE::CellType::quad4>(
          ele1, xicele1, ele1normal);
    }
    break;
    default:
      FOUR_C_THROW("rough_check_ref_config called for unknown element type");
      break;
  }
  // get center of master element
  const DRT::Element* ele2 = discret().gElement(ele2gid);
  const CORE::FE::CellType dtele2 = ele2->Shape();
  switch (dtele2)
  {
    case CORE::FE::CellType::tri3:
    {
      xicele2.PutScalar(1.0 / 3.0);
      DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::tri3>(
          ele2, xicele2, ele2coords);
    }
    break;
    case CORE::FE::CellType::tri6:
    {
      xicele2.PutScalar(1.0 / 3.0);
      DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::tri6>(
          ele2, xicele2, ele2coords);
    }
    break;
    case CORE::FE::CellType::quad4:
    {
      xicele2.PutScalar(0.0);
      DRT::UTILS::LocalToGlobalPositionAtXiRefConfig<3, CORE::FE::CellType::quad4>(
          ele2, xicele2, ele2coords);
    }
    break;
    default:
      FOUR_C_THROW("rough_check_ref_config called for unknown element type");
      break;
  }

  // calculate vector ele1 center -> ele2 center
  ele1ele2vec.Update(-1.0, ele1coords, 1.0, ele2coords, 0.0);
  // as we are only interested in the angle between connecting vector and normal we normalize it to
  // 1
  const double ele1ele2vecnorm = ele1ele2vec.Norm2();
  if (ele1ele2vecnorm > 1.0e-14) ele1ele2vec.Scale(1.0 / ele1ele2vecnorm);

  // calculate scalar product of slavemastervector and slavenormal in reference coordinates
  scalarprod.MultiplyTN(1.0, ele1ele2vec, ele1normal, 0.0);
  const double scalarprodnorm = scalarprod.MaxValue();

  // if normal (in ref. config) of ele1 and vector (in ref. config) connecting centers of ele1 and
  // ele2 point in the same direction -> integrate!
  if (scalarprodnorm > 0.0) refconfig = true;

  return refconfig;
}

/*----------------------------------------------------------------------*
 | Update and contact search (private)                     schmidt 01/19|
 *----------------------------------------------------------------------*/
void CONTACT::UnbiasedSelfBinaryTree::search_contact()
{
  // check is root node available
  if (roots().size() == 0) FOUR_C_THROW("No root node for search!");
  if (roots()[0] == Teuchos::null) FOUR_C_THROW("No root node for search!");

  // reset contact pairs from last iteration
  set_contact_pairs().clear();

  //**********************************************************************
  // STEP 1: update geometry (DOPs and sample vectors) bottom-up
  //**********************************************************************
  // update tree bottom up (for every tree layer)
  for (int i = ((int)(treenodes().size() - 1)); i >= 0; --i)
  {
    for (int j = 0; j < (int)(treenodes()[i].size()); ++j)
      treenodes()[i][j]->UpdateSlabsBottomUp(enlarge());
  }
  update_normals();

  //**********************************************************************
  // STEP 2: distribute roots among all processors
  //**********************************************************************
  // therefore find out which processor owns which root node
  std::vector<unsigned> myroots(0);
  for (unsigned i = 0; i < roots().size(); ++i)
    if (roots()[i]->Owner() == comm().MyPID()) myroots.push_back(i);

  //**********************************************************************
  // STEP 3: search for self contact starting at root nodes
  //**********************************************************************
  for (unsigned k = 0; k < myroots.size(); ++k) search_self_contact(roots()[myroots[k]]);

  //**********************************************************************
  // STEP 4: search for two-body contact between different roots
  //**********************************************************************
  // keep track of already checked root nodes
  std::set<unsigned> checkedroots;
  for (unsigned m = 0; m < myroots.size(); ++m)
  {
    checkedroots.insert(myroots[m]);
    for (unsigned k = 0; k < roots().size(); ++k)
    {
      // only perform the search if current proc did not search yet
      if (checkedroots.find(k) == checkedroots.end())
        search_root_contact(roots()[myroots[m]], roots()[k]);
    }
  }

  //**********************************************************************
  // STEP 5: all contact elements have to be slave elements
  //**********************************************************************
  std::map<int, Teuchos::RCP<SelfBinaryTreeNode>>::iterator leafiter = set_leafsmap().begin();
  std::map<int, Teuchos::RCP<SelfBinaryTreeNode>>::iterator leafiter_end = set_leafsmap().end();

  // set all contact elements and nodes to slave
  while (leafiter != leafiter_end)
  {
    const int gid = leafiter->first;
    DRT::Element* element = discret().gElement(gid);
    CONTACT::Element* celement = dynamic_cast<CONTACT::Element*>(element);

    // set contact element to slave
    celement->SetSlave() = true;

    // set nodes to slave
    for (int i = 0; i < element->num_node(); ++i)
    {
      DRT::Node* node = element->Nodes()[i];
      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
      cnode->SetSlave() = true;
    }
    // increment iterator
    ++leafiter;
  }

  //**********************************************************************
  // optional STEP 6: communicate Searchelements to all Procs
  //**********************************************************************
  if (searchele_all_proc_) communicate_search_elements_all_procs();

  // define the search elements based on the contact pairs map
  while (!contact_pairs().empty())
  {
    define_search_elements();
  }

  return;
}

/*------------------------------------------------------------------------*
 | Communicate the Search Elements to all processors (private) ager 09/19 |
 *-----------------------------------------------------------------------*/
void CONTACT::UnbiasedSelfBinaryTree::communicate_search_elements_all_procs()
{
  for (int elelid = 0; elelid < discret().ElementColMap()->NumMyElements(); ++elelid)
  {
    int elegid = discret().ElementColMap()->GID(elelid);
    std::vector<int> searchelements;
    std::vector<int> searchelements_all;
    if (contact_pairs().find(elegid) != contact_pairs().end())
      searchelements = contact_pairs()[elegid];

    CORE::LINALG::AllreduceVector(searchelements, searchelements_all, discret().Comm());

    if (searchelements_all.size()) set_contact_pairs()[elegid] = searchelements_all;
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
