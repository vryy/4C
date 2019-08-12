/*---------------------------------------------------------------------*/
/*!

\brief A collection of helper methods related to partitioning and parallel distribution

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#include "drt_utils_rebalancing.H"

#include "../linalg/linalg_utils.H"

#include <Epetra_Time.h>

/* Include Isorropia_Exception.hpp only because the helper functions at
 * the bottom of this file (which create the epetra objects) can
 * potentially throw exceptions.
 */
#include <Isorropia_Exception.hpp>

// The Isorropia symbols being demonstrated are declared in these headers:
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::REBALANCING::ComputeRebalancedNodeMaps(
    Teuchos::RCP<DRT::Discretization> discretization, Teuchos::RCP<const Epetra_Map> elementRowMap,
    Teuchos::RCP<Epetra_Map>& nodeRowMap, Teuchos::RCP<Epetra_Map>& nodeColumnMap,
    Teuchos::RCP<const Epetra_Comm> comm, const bool outflag, const int numPartitions,
    const double imbalanceTol)
{
  const int myrank = discretization->Comm().MyPID();
  const Epetra_Time timer(*comm);
  const double t1 = timer.ElapsedTime();
  if (!myrank && outflag)
    std::cout << "Rebalance nodal maps of discretization '" << discretization->Name() << "'..."
              << std::endl;

  // create nodal graph of existing problem
  Teuchos::RCP<const Epetra_CrsGraph> initialGraph =
      DRT::UTILS::REBALANCING::BuildGraph(discretization, elementRowMap, nodeRowMap, comm, outflag);

  // Create parameter list with rebalancing options
  Teuchos::RCP<Teuchos::ParameterList> rebalanceParams =
      CreateRebalancingParameterList(numPartitions, imbalanceTol);

  // Compute rebalanced graph
  Teuchos::RCP<Epetra_CrsGraph> balancedGraph =
      DRT::UTILS::REBALANCING::RebalanceGraph(initialGraph, *rebalanceParams);

  // Extract rebalanced maps
  nodeRowMap = Teuchos::rcp(new Epetra_Map(-1, balancedGraph->RowMap().NumMyElements(),
      balancedGraph->RowMap().MyGlobalElements(), 0, *comm));
  nodeColumnMap = Teuchos::rcp(new Epetra_Map(-1, balancedGraph->ColMap().NumMyElements(),
      balancedGraph->ColMap().MyGlobalElements(), 0, *comm));

  const double t2 = timer.ElapsedTime();
  if (!myrank && outflag)
    std::cout << std::setprecision(5) << "Node map rebalancing:    " << t2 - t1 << " secs"
              << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::REBALANCING::ComputeRebalancedNodeMapsUsingWeights(
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map>& rownodes,
    Teuchos::RCP<Epetra_Map>& colnodes, const bool outflag)
{
  const int myrank = dis->Comm().MyPID();
  Epetra_Time timer(dis->Comm());
  const double t1 = timer.ElapsedTime();
  if (!myrank && outflag)
    std::cout << "Rebalance nodal maps of discretization '" << dis->Name() << "'..." << std::endl;

  // create nodal graph of existing problem
  Teuchos::RCP<const Epetra_CrsGraph> initgraph = dis->BuildNodeGraph();

  // Setup cost describer based on element connectivity
  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs = SetupCostDescriber(*dis);

  // Create parameter list with repartitioning options
  Teuchos::ParameterList paramlist;
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_APPROACH", "PARTITION");

  // Compute rebalanced graph
  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
      DRT::UTILS::REBALANCING::RebalanceGraph(initgraph, costs, paramlist);

  // extract repartitioned maps
  rownodes = Teuchos::rcp(new Epetra_Map(-1, balanced_graph->RowMap().NumMyElements(),
      balanced_graph->RowMap().MyGlobalElements(), 0, dis->Comm()));
  colnodes = Teuchos::rcp(new Epetra_Map(-1, balanced_graph->ColMap().NumMyElements(),
      balanced_graph->ColMap().MyGlobalElements(), 0, dis->Comm()));

  const double t2 = timer.ElapsedTime();
  if (!myrank && outflag)
    std::cout << std::setprecision(5) << "Node map rebalancing:    " << t2 - t1 << " secs"
              << std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Isorropia::Epetra::CostDescriber> DRT::UTILS::REBALANCING::SetupCostDescriber(
    const DRT::Discretization& discretization)
{
  const Epetra_Map* oldnoderowmap = discretization.NodeRowMap();
  // Now we're going to create a Epetra_Vector with vertex weights and a Epetra_CrsMatrix
  // for the edge weights to be used in the partitioning operation.
  // weights must be at least one for zoltan
  Teuchos::RCP<Epetra_CrsMatrix> crs_ge_weights =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *oldnoderowmap, 15));
  Teuchos::RCP<Epetra_Vector> vweights = LINALG::CreateVector(*oldnoderowmap, true);

  // loop all row elements and get their cost of evaluation
  for (int i = 0; i < discretization.ElementRowMap()->NumMyElements(); ++i)
  {
    DRT::Element* ele = discretization.lRowElement(i);
    DRT::Node** nodes = ele->Nodes();
    const int numnode = ele->NumNode();
    std::vector<int> lm(numnode);
    std::vector<int> lmrowowner(numnode);
    for (int n = 0; n < numnode; ++n)
    {
      lm[n] = nodes[n]->Id();
      lmrowowner[n] = nodes[n]->Owner();
    }

    // element vector and matrix for weights of nodes and edges
    Epetra_SerialDenseMatrix edgeweigths_ele;
    Epetra_SerialDenseVector nodeweights_ele;
    // evaluate elements to get their evaluation cost
    ele->NodalConnectivity(edgeweigths_ele, nodeweights_ele);

    LINALG::Assemble(*crs_ge_weights, edgeweigths_ele, lm, lmrowowner, lm);
    LINALG::Assemble(*vweights, nodeweights_ele, lm, lmrowowner);
  }

  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs =
      Teuchos::rcp(new Isorropia::Epetra::CostDescriber);
  costs->setGraphEdgeWeights(crs_ge_weights);
  costs->setVertexWeights(vweights);

  return costs;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_CrsGraph> DRT::UTILS::REBALANCING::BuildGraph(
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<const Epetra_Map> roweles,
    Teuchos::RCP<Epetra_Map>& rownodes, Teuchos::RCP<const Epetra_Comm> comm, bool outflag)
{
  const int myrank = comm->MyPID();
  const int numproc = comm->NumProc();

  // create a set of all nodes that I have
  std::set<int> mynodes;
  for (int lid = 0; lid < roweles->NumMyElements(); ++lid)
  {
    DRT::Element* ele = dis->gElement(roweles->GID(lid));
    const int numnode = ele->NumNode();
    const int* nodeids = ele->NodeIds();
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
    comm->Broadcast(&size, 1, proc);
    if (proc != myrank) recvnodes.resize(size);
    comm->Broadcast(&recvnodes[0], size, proc);
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
    comm->Barrier();
  }


  // copy the set to a vector
  {
    std::vector<int> nodes;
    std::set<int>::iterator fool;
    for (fool = mynodes.begin(); fool != mynodes.end(); ++fool) nodes.push_back(*fool);
    mynodes.clear();
    // create a non-overlapping row map
    rownodes = Teuchos::rcp(new Epetra_Map(-1, (int)nodes.size(), &nodes[0], 0, *comm));
  }


  // start building the graph object
  std::map<int, std::set<int>> locals;
  std::map<int, std::set<int>> remotes;
  for (int lid = 0; lid < roweles->NumMyElements(); ++lid)
  {
    DRT::Element* ele = dis->gElement(roweles->GID(lid));
    const int numnode = ele->NumNode();
    const int* nodeids = ele->NodeIds();
    for (int i = 0; i < numnode; ++i)
    {
      const int lid = rownodes->LID(nodeids[i]);  // am I owner of this gid?
      std::map<int, std::set<int>>* insertmap = NULL;
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
    }  // for (int i=0; i<numnode; ++i)
  }    // for (int lid=0;lid<roweles->NumMyElements();++lid)
  // run through locals and remotes to find the max bandwith
  int maxband = 0;
  {
    int smaxband = 0;
    std::map<int, std::set<int>>::iterator fool;
    for (fool = locals.begin(); fool != locals.end(); ++fool)
      if (smaxband < (int)fool->second.size()) smaxband = (int)fool->second.size();
    for (fool = remotes.begin(); fool != remotes.end(); ++fool)
      if (smaxband < (int)fool->second.size()) smaxband = (int)fool->second.size();
    comm->MaxAll(&smaxband, &maxband, 1);
  }
  if (!myrank && outflag)
  {
    printf("parmetis max nodal bandwith %d\n", maxband);
    fflush(stdout);
  }

  Teuchos::RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *rownodes, maxband, false));
  comm->Barrier();

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
        dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d", err, grid);
    }
    locals.clear();
  }

  comm->Barrier();


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
    comm->Broadcast(&size, 1, proc);
    if (proc != myrank) recvnodes.resize(size);
    comm->Broadcast(&recvnodes[0], size, proc);
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
          if (err < 0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d", err);
          ptr += (num + 1);
        }
        else  // I don't have it so I don't care for entries of this row, goto next row
          ptr += (num + 1);
      }
    }
    comm->Barrier();
  }  //  for (int proc=0; proc<numproc; ++proc)
  remotes.clear();

  comm->Barrier();

  // finish graph
  graph->FillComplete();
  graph->OptimizeStorage();

  comm->Barrier();

  return graph;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::REBALANCING::ExportAndFillCompleteDiscretization(
    DRT::Discretization& discretization, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
    const bool assigndegreesoffreedom, const bool initelements, const bool doboundaryconditions)
{
  // Export nodes
  discretization.ExportRowNodes(noderowmap);
  discretization.ExportColumnNodes(nodecolmap);

  // Build reasonable maps for elements from the already valid and final node maps
  Teuchos::RCP<Epetra_Map> elerowmap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> elecolmap = Teuchos::null;
  discretization.BuildElementRowColumn(noderowmap, nodecolmap, elerowmap, elecolmap);
  discretization.ExportRowElements(*elerowmap);
  discretization.ExportColumnElements(*elecolmap);

  discretization.FillComplete(assigndegreesoffreedom, initelements, doboundaryconditions);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::REBALANCING::RedistributeAndFillCompleteDiscretizationUsingWeights(
    Teuchos::RCP<DRT::Discretization> discretization, const bool assigndegreesoffreedom,
    const bool initelements, const bool doboundaryconditions)
{
  // maps to be filled with final distributed node maps
  Teuchos::RCP<Epetra_Map> rownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> colnodes = Teuchos::null;

  // do weighted repartitioning to obtain new row/column maps
  DRT::UTILS::REBALANCING::ComputeRebalancedNodeMapsUsingWeights(
      discretization, rownodes, colnodes, true);

  // rebuild the discretization with new maps
  DRT::UTILS::REBALANCING::ExportAndFillCompleteDiscretization(*discretization, *rownodes,
      *colnodes, assigndegreesoffreedom, initelements, doboundaryconditions);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> DRT::UTILS::REBALANCING::RebalanceGraph(
    Teuchos::RCP<const Epetra_CrsGraph> initialGraph, const Teuchos::ParameterList& rebalanceParams)
{
  Epetra_CrsGraph* balancedGraph = NULL;
  try
  {
    balancedGraph = Isorropia::Epetra::createBalancedCopy(*initialGraph, rebalanceParams);
  }
  catch (std::exception& exc)
  {
    std::cout << "Isorropia::createBalancedCopy threw "
              << "exception '" << exc.what() << "' on proc " << initialGraph->Comm().MyPID()
              << std::endl;
    dserror("Error within Isorropia (graph balancing)");
  }

  balancedGraph->FillComplete();
  balancedGraph->OptimizeStorage();

  return Teuchos::rcp(balancedGraph);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> DRT::UTILS::REBALANCING::RebalanceGraph(
    Teuchos::RCP<const Epetra_CrsGraph> initialGraph,
    Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs,
    const Teuchos::ParameterList& rebalanceParams)
{
  Teuchos::RCP<Epetra_CrsGraph> balanced_graph = Teuchos::null;
  try
  {
    Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
        Teuchos::rcp(new Isorropia::Epetra::Partitioner(initialGraph, costs, rebalanceParams));

    Isorropia::Epetra::Redistributor rd(partitioner);
    balanced_graph = rd.redistribute(*initialGraph, true);
  }
  catch (std::exception& exc)
  {
    std::cout << "Isorropia threw exception '" << exc.what() << "' on proc "
              << initialGraph->Comm().MyPID() << std::endl;
    dserror("Error within Isorropia (graph balancing)");
  }

  return balanced_graph;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> DRT::UTILS::REBALANCING::CreateRebalancingParameterList(
    const int numPartitions, const double imbalanceTol)
{
  Teuchos::RCP<Teuchos::ParameterList> rebalancingParams =
      Teuchos::rcp(new Teuchos::ParameterList());

  rebalancingParams->set<std::string>("num parts", std::to_string(numPartitions));
  rebalancingParams->set<std::string>("IMBALANCE_TOL", std::to_string(imbalanceTol));

  return rebalancingParams;
}
