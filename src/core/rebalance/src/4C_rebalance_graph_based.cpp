/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of functions related to partitioning and parallel distribution

\level 0

*/
/*----------------------------------------------------------------------*/

#include "4C_rebalance_graph_based.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_geometric_search_bounding_volume.hpp"
#include "4C_fem_geometric_search_distributed_tree.hpp"
#include "4C_fem_geometric_search_params.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <Epetra_FECrsGraph.h>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_Exception.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::pair<Teuchos::RCP<Epetra_Map>, Teuchos::RCP<Epetra_Map>> Core::Rebalance::RebalanceNodeMaps(
    Teuchos::RCP<const Epetra_CrsGraph> initialGraph, const Teuchos::ParameterList& rebalanceParams,
    const Teuchos::RCP<Epetra_Vector>& initialNodeWeights,
    const Teuchos::RCP<Epetra_CrsMatrix>& initialEdgeWeights,
    const Teuchos::RCP<Epetra_MultiVector>& initialNodeCoordinates)
{
  TEUCHOS_FUNC_TIME_MONITOR("Rebalance::RebalanceNodeMaps");

  // Compute rebalanced graph
  Teuchos::RCP<Epetra_CrsGraph> balanced_graph = Rebalance::RebalanceGraph(*initialGraph,
      rebalanceParams, initialNodeWeights, initialEdgeWeights, initialNodeCoordinates);

  // extract repartitioned maps
  Teuchos::RCP<Epetra_Map> rownodes =
      Teuchos::rcp(new Epetra_Map(-1, balanced_graph->RowMap().NumMyElements(),
          balanced_graph->RowMap().MyGlobalElements(), 0, initialGraph->Comm()));
  Teuchos::RCP<Epetra_Map> colnodes =
      Teuchos::rcp(new Epetra_Map(-1, balanced_graph->ColMap().NumMyElements(),
          balanced_graph->ColMap().MyGlobalElements(), 0, initialGraph->Comm()));

  return {rownodes, colnodes};
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> Core::Rebalance::RebalanceGraph(const Epetra_CrsGraph& initialGraph,
    const Teuchos::ParameterList& rebalanceParams,
    const Teuchos::RCP<Epetra_Vector>& initialNodeWeights,
    const Teuchos::RCP<Epetra_CrsMatrix>& initialEdgeWeights,
    const Teuchos::RCP<Epetra_MultiVector>& initialNodeCoordinates)
{
  TEUCHOS_FUNC_TIME_MONITOR("Rebalance::RebalanceGraph");

  Isorropia::Epetra::CostDescriber costs = Isorropia::Epetra::CostDescriber();
  if (initialNodeWeights != Teuchos::null) costs.setVertexWeights(initialNodeWeights);
  if (initialEdgeWeights != Teuchos::null) costs.setGraphEdgeWeights(initialEdgeWeights);

  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner;
  if (initialNodeCoordinates != Teuchos::null)
  {
    partitioner = Teuchos::rcp(new Isorropia::Epetra::Partitioner(
        &initialGraph, &costs, initialNodeCoordinates.get(), nullptr, rebalanceParams));
  }
  else
  {
    partitioner =
        Teuchos::rcp(new Isorropia::Epetra::Partitioner(&initialGraph, &costs, rebalanceParams));
  }

  Isorropia::Epetra::Redistributor rd(partitioner);
  Teuchos::RCP<Epetra_CrsGraph> balancedGraph = rd.redistribute(initialGraph, true);

  balancedGraph->FillComplete();
  balancedGraph->OptimizeStorage();

  return balancedGraph;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::pair<Teuchos::RCP<Epetra_MultiVector>, Teuchos::RCP<Epetra_MultiVector>>
Core::Rebalance::RebalanceCoordinates(const Epetra_MultiVector& initialCoordinates,
    const Teuchos::ParameterList& rebalanceParams, const Epetra_MultiVector& initialWeights)
{
  TEUCHOS_FUNC_TIME_MONITOR("Rebalance::RebalanceCoordinates");

  Teuchos::RCP<Isorropia::Epetra::Partitioner> part = Teuchos::rcp(
      new Isorropia::Epetra::Partitioner(&initialCoordinates, &initialWeights, rebalanceParams));

  Isorropia::Epetra::Redistributor rd(part);

  return {rd.redistribute(initialCoordinates), rd.redistribute(initialWeights)};
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::pair<Teuchos::RCP<Epetra_Vector>, Teuchos::RCP<Epetra_CrsMatrix>>
Core::Rebalance::BuildWeights(const Core::FE::Discretization& dis)
{
  const Epetra_Map* noderowmap = dis.node_row_map();

  Teuchos::RCP<Epetra_CrsMatrix> crs_ge_weights =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *noderowmap, 15));
  Teuchos::RCP<Epetra_Vector> vweights = Core::LinAlg::CreateVector(*noderowmap, true);

  // loop all row elements and get their cost of evaluation
  for (int i = 0; i < dis.element_row_map()->NumMyElements(); ++i)
  {
    Core::Elements::Element* ele = dis.l_row_element(i);
    Core::Nodes::Node** nodes = ele->nodes();
    const int numnode = ele->num_node();
    std::vector<int> lm(numnode);
    std::vector<int> lmrowowner(numnode);
    for (int n = 0; n < numnode; ++n)
    {
      lm[n] = nodes[n]->id();
      lmrowowner[n] = nodes[n]->owner();
    }

    // element vector and matrix for weights of nodes and edges
    Core::LinAlg::SerialDenseMatrix edgeweigths_ele;
    Core::LinAlg::SerialDenseVector nodeweights_ele;

    // evaluate elements to get their evaluation cost
    ele->nodal_connectivity(edgeweigths_ele, nodeweights_ele);

    Core::LinAlg::Assemble(*crs_ge_weights, edgeweigths_ele, lm, lmrowowner, lm);
    Core::LinAlg::Assemble(*vweights, nodeweights_ele, lm, lmrowowner);
  }

  return {vweights, crs_ge_weights};
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_CrsGraph> Core::Rebalance::BuildGraph(
    Teuchos::RCP<Core::FE::Discretization> dis, Teuchos::RCP<const Epetra_Map> roweles)
{
  const int myrank = dis->get_comm().MyPID();
  const int numproc = dis->get_comm().NumProc();

  // create a set of all nodes that I have
  std::set<int> mynodes;
  for (int lid = 0; lid < roweles->NumMyElements(); ++lid)
  {
    Core::Elements::Element* ele = dis->g_element(roweles->GID(lid));
    const int numnode = ele->num_node();
    const int* nodeids = ele->node_ids();
    copy(nodeids, nodeids + numnode, inserter(mynodes, mynodes.begin()));
  }

  // build a unique row map from the overlapping sets
  for (int proc = 0; proc < numproc; ++proc)
  {
    int size = 0;
    std::vector<int> recvnodes;
    if (proc == myrank)
    {
      recvnodes.clear();
      std::set<int>::iterator fool;
      for (fool = mynodes.begin(); fool != mynodes.end(); ++fool) recvnodes.push_back(*fool);
      size = (int)recvnodes.size();
    }
    dis->get_comm().Broadcast(&size, 1, proc);
    if (proc != myrank) recvnodes.resize(size);
    dis->get_comm().Broadcast(&recvnodes[0], size, proc);
    if (proc != myrank)
    {
      for (int i = 0; i < size; ++i)
      {
        std::set<int>::iterator fool = mynodes.find(recvnodes[i]);
        if (fool == mynodes.end())
          continue;
        else
          mynodes.erase(fool);
      }
    }
    dis->get_comm().Barrier();
  }

  Teuchos::RCP<Epetra_Map> rownodes = Teuchos::null;
  // copy the set to a vector
  {
    std::vector<int> nodes;
    std::set<int>::iterator fool;
    for (fool = mynodes.begin(); fool != mynodes.end(); ++fool) nodes.push_back(*fool);
    mynodes.clear();
    // create a non-overlapping row map
    rownodes = Teuchos::rcp(new Epetra_Map(-1, (int)nodes.size(), &nodes[0], 0, dis->get_comm()));
  }

  // start building the graph object
  std::map<int, std::set<int>> locals;
  std::map<int, std::set<int>> remotes;
  for (int lid = 0; lid < roweles->NumMyElements(); ++lid)
  {
    Core::Elements::Element* ele = dis->g_element(roweles->GID(lid));
    const int numnode = ele->num_node();
    const int* nodeids = ele->node_ids();
    for (int i = 0; i < numnode; ++i)
    {
      const int lid = rownodes->LID(nodeids[i]);  // am I owner of this gid?
      std::map<int, std::set<int>>* insertmap = nullptr;
      if (lid != -1)
        insertmap = &locals;
      else
        insertmap = &remotes;
      // see whether we already have an entry for nodeids[i]
      std::map<int, std::set<int>>::iterator fool = (*insertmap).find(nodeids[i]);
      if (fool == (*insertmap).end())  // no entry in that row yet
      {
        std::set<int> tmp;
        copy(nodeids, nodeids + numnode, inserter(tmp, tmp.begin()));
        (*insertmap)[nodeids[i]] = tmp;
      }
      else
      {
        std::set<int>& imap = fool->second;
        copy(nodeids, nodeids + numnode, inserter(imap, imap.begin()));
      }
    }
  }

  // run through locals and remotes to find the max bandwith
  int maxband = 0;
  {
    int smaxband = 0;
    std::map<int, std::set<int>>::iterator fool;
    for (fool = locals.begin(); fool != locals.end(); ++fool)
      if (smaxband < (int)fool->second.size()) smaxband = (int)fool->second.size();
    for (fool = remotes.begin(); fool != remotes.end(); ++fool)
      if (smaxband < (int)fool->second.size()) smaxband = (int)fool->second.size();
    dis->get_comm().MaxAll(&smaxband, &maxband, 1);
  }

  Teuchos::RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *rownodes, maxband, false));
  dis->get_comm().Barrier();

  // fill all local entries into the graph
  {
    std::map<int, std::set<int>>::iterator fool = locals.begin();
    for (; fool != locals.end(); ++fool)
    {
      const int grid = fool->first;
      std::vector<int> cols(0, 0);
      std::set<int>::iterator setfool = fool->second.begin();
      for (; setfool != fool->second.end(); ++setfool) cols.push_back(*setfool);
      int err = graph->InsertGlobalIndices(grid, (int)cols.size(), &cols[0]);
      if (err < 0)
        FOUR_C_THROW(
            "Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d", err, grid);
    }
    locals.clear();
  }

  dis->get_comm().Barrier();

  // now we need to communicate and add the remote entries
  for (int proc = numproc - 1; proc >= 0; --proc)
  {
    std::vector<int> recvnodes;
    int size = 0;
    if (proc == myrank)
    {
      recvnodes.clear();
      std::map<int, std::set<int>>::iterator mapfool = remotes.begin();
      for (; mapfool != remotes.end(); ++mapfool)
      {
        recvnodes.push_back((int)mapfool->second.size() + 1);  // length of this entry
        recvnodes.push_back(mapfool->first);                   // global row id
        std::set<int>::iterator fool = mapfool->second.begin();
        for (; fool != mapfool->second.end(); ++fool)  // global col ids
          recvnodes.push_back(*fool);
      }
      size = (int)recvnodes.size();
    }
    dis->get_comm().Broadcast(&size, 1, proc);
    if (proc != myrank) recvnodes.resize(size);
    dis->get_comm().Broadcast(&recvnodes[0], size, proc);
    if (proc != myrank && size)
    {
      int* ptr = &recvnodes[0];
      while (ptr < &recvnodes[size - 1])
      {
        int num = *ptr;
        int grid = *(ptr + 1);
        // see whether I have grid in my row map
        if (rownodes->LID(grid) != -1)  // I have it, put stuff in my graph
        {
          int err = graph->InsertGlobalIndices(grid, num - 1, (ptr + 2));
          if (err < 0) FOUR_C_THROW("Epetra_CrsGraph::InsertGlobalIndices returned %d", err);
          ptr += (num + 1);
        }
        else  // I don't have it so I don't care for entries of this row, goto next row
          ptr += (num + 1);
      }
    }
    dis->get_comm().Barrier();
  }
  remotes.clear();

  dis->get_comm().Barrier();

  // finish graph
  graph->FillComplete();
  graph->OptimizeStorage();

  dis->get_comm().Barrier();

  return graph;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_CrsGraph> Core::Rebalance::BuildMonolithicNodeGraph(
    const Core::FE::Discretization& dis, const Core::GeometricSearch::GeometricSearchParams& params)
{
  // 1. Do a global geometric search
  Epetra_Vector zero_vector = Epetra_Vector(*(dis.dof_col_map()), true);

  std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> bounding_boxes;
  for (const auto* element : dis.my_row_element_range())
  {
    bounding_boxes.emplace_back(
        std::make_pair(element->id(), element->get_bounding_volume(dis, zero_vector, params)));
  }

  auto result = Core::GeometricSearch::GlobalCollisionSearch(
      bounding_boxes, bounding_boxes, dis.get_comm(), params.verbosity_);

  // 2. Set up a multivector which will be populated with all ghosting information,
  // i.e., the nodal connectivity of each element that collides with an element on this rank
  const int n_nodes_per_element_max = 27;  // element with highest node count is hex27
  Epetra_MultiVector node_information(*dis.element_row_map(), n_nodes_per_element_max, true);

  for (int rowele_i = 0; rowele_i < dis.num_my_row_elements(); ++rowele_i)
  {
    const auto* element = dis.l_row_element(rowele_i);
    for (int i_node = 0; i_node < element->num_node(); ++i_node)
    {
      const auto* node = element->nodes()[i_node];
      node_information.SumIntoMyValue(rowele_i, i_node, node->id());
    }
    node_information.SumIntoMyValue(rowele_i, element->num_node(), -1);
  }

  // 3. Get the connectivity information
  std::set<int> my_colliding_primitives;
  for (const auto& item : result)
  {
    my_colliding_primitives.insert(std::get<3>(item));
  }
  std::vector<int> my_colliding_primitives_vec;
  for (const auto& item : my_colliding_primitives)
  {
    my_colliding_primitives_vec.emplace_back(item);
  }
  Epetra_Map my_colliding_primitives_map(-1, my_colliding_primitives_vec.size(),
      my_colliding_primitives_vec.data(), 0, dis.get_comm());
  Epetra_MultiVector my_colliding_primitives_node_ids(
      my_colliding_primitives_map, n_nodes_per_element_max, false);
  Core::LinAlg::export_to(node_information, my_colliding_primitives_node_ids);

  // 4. Build and fill the graph with element internal connectivities
  auto my_graph = Teuchos::rcp(new Epetra_FECrsGraph(Copy, *(dis.node_row_map()), 40, false));

  for (const auto* element : dis.my_row_element_range())
  {
    for (int i_node = 0; i_node < element->num_node(); ++i_node)
    {
      const auto* node_main = element->nodes()[i_node];
      int index_main = node_main->id();
      for (int j_node = 0; j_node < element->num_node(); ++j_node)
      {
        const auto* node_inner = element->nodes()[j_node];
        int index = node_inner->id();

        int err = my_graph->InsertGlobalIndices(1, &index_main, 1, &index);
        if (err < 0)
          FOUR_C_THROW("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d", err,
              node_main->id());
      }
    }
  }

  // 5. Fill the graph with the geometric close entries
  for (const auto& [predicate_lid, predicate_gid, primitive_lid, primitive_gid, primitive_proc] :
      result)
  {
    int predicate_lid_discretization = dis.element_row_map()->LID(predicate_gid);
    if (predicate_lid_discretization < 0)
      FOUR_C_THROW("Could not find lid for predicate with gid %d on rank %d", predicate_gid,
          dis.get_comm().MyPID());
    if (predicate_lid != predicate_lid_discretization)
      FOUR_C_THROW("The ids dont match from arborx and the discretization");
    const auto* predicate = dis.g_element(predicate_gid);

    int primitive_lid_in_map = my_colliding_primitives_map.LID(primitive_gid);
    if (primitive_lid_in_map < 0) FOUR_C_THROW("Could not find lid for gid %d", primitive_gid);

    for (int i_node = 0; i_node < predicate->num_node(); ++i_node)
    {
      const auto* node_main = predicate->nodes()[i_node];
      int index_main = node_main->id();
      for (int j_node = 0; j_node < n_nodes_per_element_max; ++j_node)
      {
        // Get indices for primitive nodes
        int primitive_node_index =
            (int)(my_colliding_primitives_node_ids.Pointers()[j_node][primitive_lid_in_map]);

        if (primitive_node_index == -1)
          break;
        else
        {
          int err = my_graph->InsertGlobalIndices(1, &index_main, 1, &primitive_node_index);
          if (err < 0)
            FOUR_C_THROW("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d", err,
                node_main->id());
        }
      }
    }
  }

  my_graph->GlobalAssemble(true);
  my_graph->OptimizeStorage();

  return my_graph;
}
FOUR_C_NAMESPACE_CLOSE
