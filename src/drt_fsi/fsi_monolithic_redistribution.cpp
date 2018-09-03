/*----------------------------------------------------------------------------*/
/*!
\file fsi_monolithic_redistribution.cpp

\level 2

\maintainer Matthias Mayr

\brief Parallel Domain Redistribution of monolithic FSI
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
// Teuchos
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

// Zoltan / Isorropia
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

// baci
#include "fsi_monolithic.H"
#include "fsi_debugwriter.H"

#include "../drt_inpar/inpar_fsi.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_ale_fsi.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_adapter/ad_ale.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_str_fsi_timint_adaptive.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io.H"

#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_ale/ale_utils_mapextractor.H"

#include "../drt_lib/drt_node.H"

/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::RedistributeMonolithicGraph(
    const FSI_COUPLING coupling, const Epetra_Comm& comm)
{
  Epetra_Time timer(comm);

  const int myrank = comm.MyPID();

  /***********************/
  /* get interface nodes */
  /***********************/

  // initialize maps for row nodes
  std::map<int, DRT::Node*> structurenodes;
  std::map<int, DRT::Node*> fluidnodes;

  // initialize maps for column nodes
  std::map<int, DRT::Node*> structuregnodes;
  std::map<int, DRT::Node*> fluidgnodes;

  // initialize maps for elements
  std::map<int, Teuchos::RCP<DRT::Element>> structureelements;
  std::map<int, Teuchos::RCP<DRT::Element>> fluidelements;


  // access the discretizations
  Teuchos::RCP<DRT::Discretization> structuredis = DRT::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis = DRT::Problem::Instance()->GetDis("ale");

  // Fill maps based on condition for master side (masterdis != slavedis)
  DRT::UTILS::FindConditionObjects(
      *structuredis, structurenodes, structuregnodes, structureelements, "FSICoupling");

  // Fill maps based on condition for slave side (masterdis != slavedis)
  DRT::UTILS::FindConditionObjects(
      *fluiddis, fluidnodes, fluidgnodes, fluidelements, "FSICoupling");

  std::map<int, DRT::Node*>* slavenodesPtr = NULL;
  std::map<int, DRT::Node*>* mastergnodesPtr = NULL;

  // ToDo (mayr) Move this to routine in derived classes to replace the if-clause by inheritence
  if (coupling == fsi_iter_mortar_monolithicfluidsplit or
      coupling == fsi_iter_sliding_monolithicfluidsplit)
  {
    slavenodesPtr = &fluidnodes;
    mastergnodesPtr = &structuregnodes;
  }
  else if (coupling == fsi_iter_mortar_monolithicstructuresplit or
           coupling == fsi_iter_sliding_monolithicstructuresplit)
  {
    slavenodesPtr = &structurenodes;
    mastergnodesPtr = &fluidgnodes;
  }
  else
    dserror(
        "\nDomain redistribution / rebalancing only implemented for mortar "
        "coupling!");

  std::map<int, std::vector<int>>
      fluidToStructureMap;  // maps fluid nodes to opposing structure nodes
  std::map<int, std::vector<int>>
      structureToFluidMap;  // maps fluid nodes to opposing structure nodes

  /*********************************/
  /* get node mapping at interface */
  /*********************************/
  CreateInterfaceMapping(structuredis, fluiddis, slavenodesPtr, mastergnodesPtr,
      fluidToStructureMap, structureToFluidMap);


  /***************************/
  /* create monolithic graph */
  /***************************/

  // create nodal graphs and maps of structure and fluid field
  Teuchos::RCP<const Epetra_CrsGraph> structureGraph = structuredis->BuildNodeGraph();
  Teuchos::RCP<const Epetra_CrsGraph> fluidGraph = fluiddis->BuildNodeGraph();
  const Epetra_Map& structureGraph_map = (Epetra_Map&)structureGraph->RowMap();
  const Epetra_Map& fluidGraph_map = (Epetra_Map&)fluidGraph->RowMap();

  // Create monolithic map. Interface fluid nodes are omitted.

  int numGlobalStructureNodes = structureGraph_map.NumGlobalElements();

  int numMyStructureNodes = structureGraph_map.NumMyElements();
  int numMyFluidNodes = fluidGraph_map.NumMyElements();

  int prelimNumMyNodes = numMyStructureNodes + numMyFluidNodes;

  // vector of global gids of structure + fluid
  int gid[prelimNumMyNodes];
  int* gid_structure = structureGraph_map.MyGlobalElements();
  int* gid_fluid = fluidGraph_map.MyGlobalElements();

  for (int i = 0; i < numMyStructureNodes; ++i) gid[i] = gid_structure[i];
  int countMyFluidNodes = 0;
  for (int i = 0; i < numMyFluidNodes; ++i)
  {
    try
    {
      fluidToStructureMap.at(gid_fluid[i]);
    }
    catch (std::exception& exc)
    {
      gid[numMyStructureNodes + countMyFluidNodes] = gid_fluid[i];
      countMyFluidNodes++;
    }
  }

  int countGlobalFluidNodes = 0;
  comm.SumAll(&countMyFluidNodes, &countGlobalFluidNodes, 1);
  int numMyNodes = numMyStructureNodes + countMyFluidNodes;
  int numGlobalNodes = numGlobalStructureNodes + countGlobalFluidNodes;

  const Epetra_Map& monolithicMap = Epetra_Map(numGlobalNodes, numMyNodes, gid, 0, comm);

  // create monolithic graph
  Teuchos::RCP<Epetra_CrsGraph> monolithicGraph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, monolithicMap, 0));

  std::map<int, std::vector<int>> deletedEdges;  // edges deleted during construction of monolithic
                                                 // graph: have to be inserted after redistribution
  std::map<int, std::vector<int>>
      insertedEdges;  // edges inserted during construction of monolithic graph:
                      // have to be removed after redistribution

  BuildMonolithicGraph(monolithicGraph, deletedEdges, insertedEdges, fluidToStructureMap,
      structureToFluidMap, structuredis, fluiddis);

  if (myrank == 0) std::cout << "About to call Zoltan..." << std::endl;

  /******************/
  /* redistribution */
  /******************/
  // dummy matrix for Zoltan call function
  Teuchos::RCP<Epetra_CrsMatrix> dummy;
  Teuchos::RCP<Epetra_Vector> crs_hge_weights = Teuchos::rcp(new Epetra_Vector(monolithicMap));
  crs_hge_weights->PutScalar(1.0);
  //  Teuchos::RCP<Epetra_CrsGraph> bal_graph = CallZoltan(monolithicGraph, dummy,
  //  crs_hge_weights,"HYPERGRAPH", 0);

  Teuchos::RCP<Epetra_CrsGraph> bal_graph =
      CallZoltanWithoutWeights(monolithicGraph, "HYPERGRAPH", 0);

  if (myrank == 0) std::cout << "... returned from Zoltan." << std::endl;

  // get maps of monolithic graph with deleted fluid interface nodes and inserted couplings

  bal_graph->FillComplete();

  // Extract row map
  Teuchos::RCP<Epetra_Map> monolithicRownodes;
  monolithicRownodes = Teuchos::rcp(new Epetra_Map(
      -1, bal_graph->RowMap().NumMyElements(), bal_graph->RowMap().MyGlobalElements(), 0, comm));

  // Extract col map
  Teuchos::RCP<Epetra_Map> monolithicColnodes;
  monolithicColnodes = Teuchos::rcp(new Epetra_Map(
      -1, bal_graph->ColMap().NumMyElements(), bal_graph->ColMap().MyGlobalElements(), 0, comm));

  /*************************************************/
  /* reconstruct actual fluid and structure graphs */
  /*************************************************/

  // get fluid and structure row maps with inserted fluid interface nodes (final rowmaps)

  Teuchos::RCP<Epetra_Map> structureRowmap =
      GetRedistRowMap(structureGraph_map, monolithicRownodes, fluidToStructureMap);

  Teuchos::RCP<Epetra_Map> fluidRowmap =
      GetRedistRowMap(fluidGraph_map, monolithicRownodes, fluidToStructureMap, true);

  // Now create structure and fluid graph with inserted / removed edges such that the final column
  // maps can be extracted
  Teuchos::RCP<Epetra_CrsGraph> structureGraphRedist =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *structureRowmap, 0));
  Teuchos::RCP<Epetra_CrsGraph> fluidGraphRedist =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *fluidRowmap, 0));

  RestoreRedistStructFluidGraph(insertedEdges, deletedEdges, bal_graph, monolithicRownodes,
      monolithicColnodes, structureGraphRedist, fluidGraphRedist, fluidToStructureMap);

  structureGraphRedist->FillComplete();
  fluidGraphRedist->FillComplete();

  // Extract structure and fluid column maps
  Teuchos::RCP<Epetra_Map> structureColmap;
  structureColmap = Teuchos::rcp(new Epetra_Map(-1, structureGraphRedist->ColMap().NumMyElements(),
      structureGraphRedist->ColMap().MyGlobalElements(), 0, comm));

  Teuchos::RCP<Epetra_Map> fluidColmap;
  fluidColmap = Teuchos::rcp(new Epetra_Map(-1, fluidGraphRedist->ColMap().NumMyElements(),
      fluidGraphRedist->ColMap().MyGlobalElements(), 0, comm));

  /*************************/
  /* Actual redistribution */
  /*************************/
  /*
   * Now the actual redistribution, i.e. the redistribution of the
   * discretization is done. Some vectors and the time integrators have to
   * be rebuilt as well in order to be conforming to the new row maps.
   */

  // redistribute nodes to procs
  structuredis->Redistribute(*structureRowmap, *structureColmap);
  fluiddis->Redistribute(*fluidRowmap, *fluidColmap);
  aledis->Redistribute(*fluidRowmap, *fluidColmap);

  // Distribute dofs in time integration vectors to right procs by creating time integrators new.
  // The Control File has to be rewritten as well.
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Teuchos::ParameterList> ioflags =
      Teuchos::rcp(new Teuchos::ParameterList(problem->IOParams()));
  Teuchos::RCP<IO::DiscretizationWriter> structureoutput = structuredis->Writer();
  Teuchos::RCP<IO::DiscretizationWriter> fluidoutput = fluiddis->Writer();
  Teuchos::RCP<IO::DiscretizationWriter> aleoutput = aledis->Writer();

  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();

  fluidoutput->OverwriteResultFile();
  aleoutput->OverwriteResultFile();
  structureoutput->OverwriteResultFile();
  CreateStructureTimeIntegrator(fsidyn, structuredis);
  SetLambda();

  fluidoutput->OverwriteResultFile();
  aleoutput->OverwriteResultFile();
  structureoutput->OverwriteResultFile();
  CreateFluidAndALETimeIntegrator(fsidyn, fluiddis, aledis);
  SetLambda();


  /*
   * In the control file the three fields have to appear in the order
   * structure - fluid - ale. The easiest and most stable way to achieve this
   * is to overwrite the file and write the three fields in this order.
   */

  structureoutput->OverwriteResultFile();
  fluidoutput->OverwriteResultFile();
  aleoutput->OverwriteResultFile();

  structureoutput->WriteMesh(0, 0.0);
  fluidoutput->WriteMesh(0, 0.0);
  aleoutput->WriteMesh(0, 0.0);

  // setup has do be done again
  SetNotSetup();

  if (myrank == 0) printf("Redistribution of domain in %f seconds.\n", timer.ElapsedTime());

  // just to be safe
  comm.Barrier();

  return;
}

/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::RedistributeDomainDecomposition(const INPAR::FSI::Redistribute domain,
    const FSI_COUPLING coupling, const double inputWeight1, const double inputWeight2,
    const Epetra_Comm& comm, int unbalance)
{
  Epetra_Time timer(comm);

  const int myrank = comm.MyPID();

  interfaceprocs_.clear();

  /***********************/
  /* get interface nodes */
  /***********************/

  // initialize maps for row nodes
  std::map<int, DRT::Node*> structurenodes;
  std::map<int, DRT::Node*> fluidnodes;

  // initialize maps for column nodes
  std::map<int, DRT::Node*> structuregnodes;
  std::map<int, DRT::Node*> fluidgnodes;

  // initialize maps for elements
  std::map<int, Teuchos::RCP<DRT::Element>> structureelements;
  std::map<int, Teuchos::RCP<DRT::Element>> fluidelements;


  // access the discretizations
  Teuchos::RCP<DRT::Discretization> structuredis = DRT::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis = DRT::Problem::Instance()->GetDis("ale");

  // Fill maps based on condition for master side (masterdis != slavedis)
  DRT::UTILS::FindConditionObjects(
      *structuredis, structurenodes, structuregnodes, structureelements, "FSICoupling");

  // Fill maps based on condition for slave side (masterdis != slavedis)
  DRT::UTILS::FindConditionObjects(
      *fluiddis, fluidnodes, fluidgnodes, fluidelements, "FSICoupling");

  std::map<int, DRT::Node*>* slavenodesPtr = NULL;
  std::map<int, DRT::Node*>* mastergnodesPtr = NULL;

  // ToDo (mayr) Move this to routine in derived classes to replace the if-clause by inheritence
  if (coupling == fsi_iter_mortar_monolithicfluidsplit or
      coupling == fsi_iter_sliding_monolithicfluidsplit)
  {
    slavenodesPtr = &fluidnodes;
    mastergnodesPtr = &structuregnodes;
  }
  else if (coupling == fsi_iter_mortar_monolithicstructuresplit or
           coupling == fsi_iter_sliding_monolithicstructuresplit)
  {
    slavenodesPtr = &structurenodes;
    mastergnodesPtr = &fluidgnodes;
  }
  else
    dserror(
        "\nDomain redistribution / rebalancing only implemented for mortar "
        "coupling!");


  /*******************************************/
  /* distribute masternodes to future owners */
  /*******************************************/

  std::map<int, int> nodeOwner;  // maps global node id and the owner of the neighbouring node
  std::map<int, int>* nodeOwnerPtr = &nodeOwner;
  std::map<int, std::list<int>> inverseNodeOwner;  // Create a map that maps Owner <-> Nodes
  std::map<int, std::list<int>>* inverseNodeOwnerPtr = &inverseNodeOwner;

  CreateNodeOwnerRelationship(nodeOwnerPtr, inverseNodeOwnerPtr, slavenodesPtr, mastergnodesPtr,
      structuredis, fluiddis, domain);

  /*********************************************/
  /* build nodal graph and insert edge weights */
  /*********************************************/

  Teuchos::RCP<DRT::Discretization> dis;
  if (domain == INPAR::FSI::Redistribute_structure)
    dis = structuredis;
  else if (domain == INPAR::FSI::Redistribute_fluid)
    dis = fluiddis;

  // create nodal graph of master field
  Teuchos::RCP<const Epetra_CrsGraph> initgraph = dis->BuildNodeGraph();
  const Epetra_Map& initgraph_map = (Epetra_Map&)initgraph->RowMap();

  Teuchos::RCP<Epetra_CrsMatrix> crs_ge_weights =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, initgraph_map, 0));
  Teuchos::RCP<Epetra_CrsGraph> initgraph_manip =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, initgraph_map, 0));


  std::map<int, std::list<int>> deletedEdges;
  std::map<int, std::list<int>>* deletedEdgesPtr = &deletedEdges;

  BuildWeightedGraph(crs_ge_weights, initgraph_manip, initgraph, inputWeight1, inputWeight2,
      nodeOwnerPtr, inverseNodeOwnerPtr, deletedEdgesPtr, comm);

  /******************/
  /* redistribution */
  /******************/

  // dummy vector for Zoltan call function
  Teuchos::RCP<Epetra_Vector> dummy;
  //  Teuchos::RCP<Epetra_CrsGraph> bal_graph = CallZoltan(initgraph_manip,
  //  crs_ge_weights,dummy,"GRAPH",unbalance);

  Teuchos::RCP<Epetra_CrsGraph> bal_graph = CallZoltanWithoutWeights(initgraph_manip, "GRAPH", 0);

  bal_graph->FillComplete();
  bal_graph->OptimizeStorage();

  // Extract row map
  Teuchos::RCP<Epetra_Map> rownodes;
  rownodes = Teuchos::rcp(new Epetra_Map(
      -1, bal_graph->RowMap().NumMyElements(), bal_graph->RowMap().MyGlobalElements(), 0, comm));

  /******************/
  /* switch domains */
  /******************/
  /*
   * The distribution of the nodes to procs is done. Now we want to
   * switch these domains of nodes among the procs (if necessary) such
   * that at the interface fluid and structure nodes are handled by the
   * same proc.
   */

  Teuchos::RCP<Epetra_CrsGraph> switched_bal_graph =
      SwitchDomains(rownodes, nodeOwnerPtr, bal_graph, comm);

  // Extract row map from the final graph. Column map will be done in a while, after we insert the
  // deleted edges from above.
  Teuchos::RCP<Epetra_Map> switched_rownodes;
  switched_rownodes = Teuchos::rcp(new Epetra_Map(-1, switched_bal_graph->RowMap().NumMyElements(),
      switched_bal_graph->RowMap().MyGlobalElements(), 0, comm));


  // TODO: get only rowmap from switchdomains and don't build graph.

  /*
   * We inserted edges between all interface nodes on one proc in BuildWeightedGraph
   * which would lead to an incorrect column map if we used switched_bal_graph
   * to extract the column map. Therefore, we export the initgraph to the final
   * switched_rownodes and get the columnmap afterwards.
   */


  Teuchos::RCP<Epetra_CrsGraph> final_graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *switched_rownodes, 0));
  Epetra_Export exporter(initgraph_map, *switched_rownodes);

  int err = final_graph->Export(*initgraph, exporter, Insert);
  if (err) dserror("Export not successful, error code %d!", err);

  InsertDeletedEdges(deletedEdgesPtr, switched_rownodes, final_graph);

  // Now extract column map from the final graph.
  Teuchos::RCP<Epetra_Map> switched_colnodes;
  switched_colnodes = Teuchos::rcp(new Epetra_Map(-1, final_graph->ColMap().NumMyElements(),
      final_graph->ColMap().MyGlobalElements(), 0, comm));


  /*************************/
  /* Actual redistribution */
  /*************************/
  /*
   * Now the actual redistribution, i.e. the redistribution of the
   * discretization is done. Some vectors and the time integrators have to
   * be rebuilt as well in order to be conforming to the new row maps.
   */

  // redistribute nodes to procs
  dis->Redistribute(*switched_rownodes, *switched_colnodes);

  if (domain == INPAR::FSI::Redistribute_fluid)
  {
    aledis->Redistribute(*switched_rownodes, *switched_colnodes);
  }

  // Distribute dofs in time integration vectors to right procs by creating time integrators new.
  // The Control File has to be rewritten as well.
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Teuchos::ParameterList> ioflags =
      Teuchos::rcp(new Teuchos::ParameterList(problem->IOParams()));
  Teuchos::RCP<IO::DiscretizationWriter> structureoutput = structuredis->Writer();
  Teuchos::RCP<IO::DiscretizationWriter> fluidoutput = fluiddis->Writer();
  Teuchos::RCP<IO::DiscretizationWriter> aleoutput = aledis->Writer();

  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();

  if (domain == INPAR::FSI::Redistribute_structure)
  {
    fluidoutput->OverwriteResultFile();
    aleoutput->OverwriteResultFile();
    structureoutput->OverwriteResultFile();
    CreateStructureTimeIntegrator(fsidyn, structuredis);
    SetLambda();
  }
  else if (domain == INPAR::FSI::Redistribute_fluid)
  {
    fluidoutput->OverwriteResultFile();
    aleoutput->OverwriteResultFile();
    structureoutput->OverwriteResultFile();
    CreateFluidAndALETimeIntegrator(fsidyn, fluiddis, aledis);
    SetLambda();
  }

  //  // add missing sections "fluid" and "ale" in control file
  //  if (DRT::INPUT::IntegralValue<int>(*ioflags,"OUTPUT_BIN")){
  //    if (domain == INPAR::FSI::Redistribute_structure){
  //      fluidoutput->WriteMesh(0, 0.0);
  //      aleoutput->WriteMesh(0, 0.0);
  //    }
  //    else if (domain == INPAR::FSI::Redistribute_fluid){
  ////      structureoutput->OverwriteResultFile();
  ////      fluidoutput->OverwriteResultFile();
  ////      aleoutput->OverwriteResultFile();
  //      structureoutput->WriteMesh(0, 0.0);
  ////      fluidoutput->WriteMesh(0, 0.0);
  ////      aleoutput->WriteMesh(0, 0.0);
  //    }
  //  }


  //  fluidoutput->OverwriteResultFile();
  //  aleoutput->OverwriteResultFile();
  //  structureoutput->OverwriteResultFile();
  //
  //  CreateStructureTimeIntegrator(fsidyn, structuredis);
  //  CreateFluidAndALETimeIntegrator(fsidyn, fluiddis, aledis);
  //  SetLambda();

  /*
   * In the control file the three fields have to appear in the order
   * structure - fluid - ale. The easiest and most stable way to achieve this
   * is to overwrite the file and write the three fields in this order.
   */

  structureoutput->OverwriteResultFile();
  fluidoutput->OverwriteResultFile();
  aleoutput->OverwriteResultFile();

  structureoutput->WriteMesh(0, 0.0);
  fluidoutput->WriteMesh(0, 0.0);
  aleoutput->WriteMesh(0, 0.0);


  /*
   * TODO: make sure that nodeowner is written only in first time step
   */


  // check and fix ml nullspace if neccessary

  //  // get parameter list of structural dynamics
  //  const Teuchos::ParameterList& sdyn = problem->StructuralDynamicParams();
  //  // get the solver number used for structural solver
  //  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  //  Teuchos::ParameterList slist = problem->SolverParams(slinsolvernumber);
  //  structuredis->ComputeNullSpaceIfNecessary(slist,true);
  //
  //  // get parameter list of fluid dynamics
  //  const Teuchos::ParameterList& fdyn = problem->FluidDynamicParams();
  //  // get the solver number used for fluid solver
  //  const int flinsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  //  Teuchos::ParameterList flist = problem->SolverParams(flinsolvernumber);
  //  fluiddis->ComputeNullSpaceIfNecessary(flist,true);
  //
  //  // get parameter list of ale dynamics
  //  const Teuchos::ParameterList& adyn = problem->AleDynamicParams();
  //  // get the solver number used for ale solver
  //  const int alinsolvernumber = adyn.get<int>("LINEAR_SOLVER");
  //  Teuchos::ParameterList alist = problem->SolverParams(alinsolvernumber);
  //  aledis->ComputeNullSpaceIfNecessary(alist,true);


  // setup has do be done again
  SetNotSetup();

  if (myrank == 0) printf("Redistribution of domain in %f seconds.\n", timer.ElapsedTime());

  comm.Barrier();

  return;
}

/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::BuildWeightedGraph(Teuchos::RCP<Epetra_CrsMatrix> crs_ge_weights,
    Teuchos::RCP<Epetra_CrsGraph> initgraph_manip, Teuchos::RCP<const Epetra_CrsGraph> initgraph,
    const double inputWeight1, const double inputWeight2, std::map<int, int>* nodeOwner,
    std::map<int, std::list<int>>* inverseNodeOwner, std::map<int, std::list<int>>* deletedEdges,
    const Epetra_Comm& comm)
{
  /*********************************************/
  /* build nodal graph and insert edge weights */
  /*********************************************/

  int numproc = comm.NumProc();
  int myrank = comm.MyPID();

  const Epetra_Map& initgraph_map = (Epetra_Map&)initgraph->RowMap();

  int owner;   // owner in first loop
  int owner2;  // owner in second loop
  int noderowcopy;
  int ownercopy;

  const double weight = inputWeight2;
  const double zero = 0.001;
  int graphLength = crs_ge_weights->NumGlobalCols();
  int graphIndices[graphLength];
  int graphIndicesCopy[graphLength];
  int numGraphEntries;
  int numGraphEntriesCopy;

  int numNodes = (int)nodeOwner->size();
  int maxNodes;
  comm.MaxAll(&numNodes, &maxNodes, 1);

  int numMyGraphRows = initgraph->NumMyRows();
  int maxNumGraphRows;
  comm.MaxAll(&numMyGraphRows, &maxNumGraphRows, 1);

  for (int i = 0; i < maxNumGraphRows; ++i)
  {  // all procs have to be in the loop until the very last checked all its nodes

    int graphRowGID = crs_ge_weights->GRID(i);
    if (graphRowGID != -1)
    {
      int success =
          initgraph->ExtractGlobalRowCopy(graphRowGID, graphLength, numGraphEntries, graphIndices);
      if (success != 0) dserror("\nExtract global values failed, error code %d!", success);
    }

    for (int proc = 0; proc < numproc; ++proc)
    {  // each processor checks if graphRowGID is an interface node

      int foundInterfaceNode = 0;
      noderowcopy = graphRowGID;
      comm.Broadcast(&noderowcopy, 1, proc);
      numGraphEntriesCopy = numGraphEntries;
      comm.Broadcast(&numGraphEntriesCopy, 1, proc);

      if (noderowcopy != -1)
      {  // only consider procs that still have graph rows in their map

        owner = -1;
        try
        {                                      // check which nodes are interface nodes
          owner = nodeOwner->at(noderowcopy);  // and if so: get (desired) owner
          foundInterfaceNode = 1;
        }
        catch (std::exception& exc)  // not an interface node
        {
        }

        int foundInterfaceNodeSum;
        comm.SumAll(&foundInterfaceNode, &foundInterfaceNodeSum, 1);

        if (myrank == proc)
        {
          for (int c = 0; c < numGraphEntries; ++c) graphIndicesCopy[c] = graphIndices[c];
        }

        comm.Broadcast(&graphIndicesCopy[0], numGraphEntriesCopy, proc);

        if (foundInterfaceNodeSum == 0)
        {  // not an interface node
          if (initgraph_map.LID(noderowcopy) != -1)
          {
            double fieldWeight[numGraphEntriesCopy];
            for (int w = 0; w < numGraphEntriesCopy; ++w)
            {
              fieldWeight[w] = inputWeight1;
            }

            int success;
            success = crs_ge_weights->InsertGlobalValues(
                noderowcopy, numGraphEntriesCopy, &fieldWeight[0], &graphIndicesCopy[0]);
            if (success != 0) dserror("\nInsert global values failed, error code %d!", success);
            success = initgraph_manip->InsertGlobalIndices(
                noderowcopy, numGraphEntriesCopy, &graphIndicesCopy[0]);
            if (success != 0) dserror("\nInsert global indices failed, error code %d!", success);
          }
        }
        else if (foundInterfaceNodeSum == 1)
        {  // interface nodes

          for (int proc2 = 0; proc2 < numproc; ++proc2)
          {
            ownercopy = owner;
            comm.Broadcast(&ownercopy, 1, proc2);

            if (ownercopy != -1)
            {
              // we want to connect all interface nodes on one future proc with high weights
              std::list<int> interfaceNodesOnProc = (*inverseNodeOwner)[ownercopy];
              int listLength = interfaceNodesOnProc.size();

              // find interface node pairs

              double GIDweight;  // weight for graph position noderowcopy - graphIndices[k]
              double insertWeights[listLength];
              int insertIndices[listLength];

              int insert_k = 0;

              for (int k = 0; k < numGraphEntriesCopy; ++k)
              {
                GIDweight = -1.0;
                int found = 0;
                owner2 = -1;

                if (graphIndicesCopy[k] != noderowcopy)
                {
                  try
                  {  // check which nodes are interface nodes
                    owner2 = nodeOwner->at(graphIndicesCopy[k]);  // and if so: get (desired) owner
                    if (owner2 == ownercopy)
                    {  // node wants to be on the same proc as node related to noderowcopy
                      GIDweight = weight;
                    }
                    else
                    {  // border between two neighboring patches
                      GIDweight = zero;
                    }
                    found = 1;
                  }
                  catch (std::exception& exc)  // not an interface node
                  {
                  }
                }

                int foundSum;
                comm.SumAll(&found, &foundSum, 1);

                if (foundSum == 1)
                {  // someone found an interface node

                  double GIDweightCopy;
                  for (int p = 0; p < numproc; ++p)
                  {  // collect weight info in insertWeights and insertIndices
                    GIDweightCopy = GIDweight;
                    comm.Broadcast(&GIDweightCopy, 1, p);
                    if (GIDweightCopy != -1)
                    {
                      if (GIDweightCopy != zero)
                      {
                        insertWeights[insert_k] = GIDweightCopy;
                        insertIndices[insert_k] = graphIndicesCopy[k];
                        interfaceNodesOnProc.remove(graphIndicesCopy[k]);
                        insert_k++;
                      }
                      if (GIDweightCopy == zero)
                      {  // delete that edge
                        // deletedEdges[noderowcopy] = graphIndicesCopy[k];
                        std::list<int>::iterator it = (*deletedEdges)[noderowcopy].begin();
                        (*deletedEdges)[noderowcopy].insert(it, graphIndicesCopy[k]);
                        interfaceNodesOnProc.remove(graphIndicesCopy[k]);
                      }
                      break;
                    }
                  }
                }
                else if (foundSum == 0)
                {  // no one found an interface node
                  insertWeights[insert_k] = inputWeight1;
                  insertIndices[insert_k] = graphIndicesCopy[k];
                  insert_k++;
                }
                else
                {  // we expect nodeOwner to be one to one, i.e. foundSum =  {0, 1}
                  dserror(
                      "\nNode <-> Owner mapping not one-to-one, %d procs have the same node in "
                      "their maps!",
                      foundSum);
                }
              }


              // connect interface nodes on one proc with high weights
              std::list<int>::iterator listIt = interfaceNodesOnProc.begin();
              listLength = interfaceNodesOnProc.size();

              for (int l = 0; l < listLength; ++l)
              {
                insertIndices[insert_k] = *listIt;
                insertWeights[insert_k] = weight;
                insert_k++;
                listIt++;
              }


              if (initgraph_map.LID(noderowcopy) != -1)
              {  // get processor which accesses row noderowcopy and insert weights
                int success;
                success = crs_ge_weights->InsertGlobalValues(
                    noderowcopy, insert_k, &insertWeights[0], &insertIndices[0]);
                if (success != 0)
                  dserror("\nReplace global values failed, error code %d!", success);
                success =
                    initgraph_manip->InsertGlobalIndices(noderowcopy, insert_k, &insertIndices[0]);
                if (success != 0)
                  dserror("\nInsert global indices failed, error code %d!", success);
              }
            }  // if (ownercopy != -1){
          }    // for (int proc2=0; proc2<numproc; ++proc2){
        }
        else
          dserror(
              "\n\nNode <-> Owner mapping not one-to-one, %d procs have the same node in their "
              "maps!",
              foundInterfaceNodeSum);

      }  // if (noderowcopy != -1){

      comm.Barrier();

    }  // for (int proc=0; proc<numproc; ++proc){

  }  // for (int i=0; i<maxNodes; ++i){

  crs_ge_weights->FillComplete();
  crs_ge_weights->OptimizeStorage();

  initgraph_manip->FillComplete();
  initgraph_manip->OptimizeStorage();
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> FSI::BlockMonolithic::CallZoltan(
    Teuchos::RCP<const Epetra_CrsGraph> initgraph_manip,
    Teuchos::RCP<const Epetra_CrsMatrix> matrix_weights, Teuchos::RCP<Epetra_Vector> vector_weights,
    std::string partitioningMethod, int unbalance)
{
  dserror("You should not end up here.");

  Teuchos::ParameterList paramlist;

  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs =
      Teuchos::rcp(new Isorropia::Epetra::CostDescriber);

  if (partitioningMethod == "GRAPH") costs->setGraphEdgeWeights(matrix_weights);

  if (partitioningMethod == "HYPERGRAPH") costs->setHypergraphEdgeWeights(vector_weights);

  int numproc = initgraph_manip->Comm().NumProc();
  std::stringstream ss;
  ss << numproc - unbalance;

  paramlist.set("PARTITIONING METHOD", partitioningMethod);
  paramlist.set("PRINT ZOLTAN METRICS", "2");
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("GRAPH_PACKAGE", "PHG");
  sublist.set("EDGE_WEIGHT_DIM", "1");  // One weight per edge
  sublist.set(
      "LB_APPROACH", "PARTITION");  // Build partition from scratch, in contrast to "REPARTITION"
  sublist.set("NUM_GLOBAL_PARTS", ss.str());

  //  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs = Teuchos::rcp(
  //      new Isorropia::Epetra::CostDescriber);
  //
  //  costs->setGraphEdgeWeights(crs_ge_weights);

  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
      Teuchos::rcp(new Isorropia::Epetra::Partitioner(initgraph_manip, costs, paramlist));

  Isorropia::Epetra::Redistributor rd(partitioner);
  Teuchos::RCP<Epetra_CrsGraph> bal_graph;

  /* Use a try-catch block because Isorropia will throw an exception
  if it encounters an error. */
  try
  {
    bal_graph = rd.redistribute(*initgraph_manip, false);
  }
  catch (std::exception& exc)
  {
    std::cout << "Redistribute domain: Isorropia::Epetra::Redistributor threw "
              << "exception '" << exc.what() << std::endl;
    MPI_Finalize();
  }

  // bal_graph->FillComplete();
  // bal_graph->OptimizeStorage();

  return bal_graph;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> FSI::BlockMonolithic::CallZoltanWithoutWeights(
    Teuchos::RCP<const Epetra_CrsGraph> initgraph_manip, std::string partitioningMethod,
    int unbalance)
{
  int numproc = initgraph_manip->Comm().NumProc();
  std::stringstream ss;
  ss << numproc - unbalance;
  const int parts = numproc - unbalance;

  //  paramlist.set("PARTITIONING METHOD", partitioningMethod);
  //  paramlist.set("PRINT ZOLTAN METRICS", "2");
  //  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  //  sublist.set("GRAPH_PACKAGE", "PHG");
  //  sublist.set("EDGE_WEIGHT_DIM", "1");      // One weight per edge
  //  sublist.set("LB_APPROACH", "PARTITION");  // Build partition from scratch, in contrast to
  //  "REPARTITION" sublist.set("NUM_GLOBAL_PARTS", ss.str());

  Teuchos::ParameterList paramlist;
  // No parameters. By default, Isorropia will use Zoltan hypergraph
  // partitioning, treating the graph columns as hyperedges and the
  // graph rows as vertices.

  // debug (mayr) Set imbalance tolerance to avoid empty procs
  paramlist.set("IMBALANCE_TOL", "1.03");
  //  paramlist.set("partitioning method", "hypergraph");
  //  paramlist.set("PARTITIONING METHOD", "HYPERGRAPH");
  //
  //  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  //  sublist.set("GRAPH_PACKAGE", "PHG");
  //  sublist.set("HYPERGRAPH_PACKAGE", "PHG");
  //  sublist.set("EDGE_WEIGHT_DIM", "1");      // One weight per edge
  //  sublist.set("LB_APPROACH", "PARTITION");  // Build partition from scratch, in contrast to
  //  "REPARTITION"
  ////  sublist.set("LB_METHOD", "HYPERGRAPH");
  //  sublist.set("NUM_GLOBAL_PARTS", ss.str());

  // if the user wants to use less procs than available (as indicated by
  // the input flag "parts" above) then pass on this information to the
  // parameter list for Zoltan/Isorropia
  if (parts != -1)
  {
    std::stringstream ss;
    ss << parts;
    std::string s = ss.str();
    paramlist.set("num parts", s);
  }

  //  Epetra_CrsGraph *bal_graph = NULL;
  //  try {
  //    bal_graph =
  //      Isorropia::Epetra::createBalancedCopy(*initgraph_manip, paramlist);
  //
  //  }
  //  catch(std::exception& exc) {
  //    std::cout << "Isorropia::createBalancedCopy threw "
  //         << "exception '" << exc.what() << "' on proc "
  //         << initgraph_manip->Comm().MyPID() << std::endl;
  //    dserror("Error within Isorropia (graph balancing)");
  //  }
  //
  //  Teuchos::RCP<Epetra_CrsGraph> rcp_balanced_graph = Teuchos::rcp(bal_graph);


  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
      Teuchos::rcp(new Isorropia::Epetra::Partitioner(initgraph_manip, paramlist, true));
  Teuchos::RCP<Isorropia::Epetra::Redistributor> rd =
      Teuchos::rcp(new Isorropia::Epetra::Redistributor(partitioner));
  Teuchos::RCP<Epetra_CrsGraph> rcp_balanced_graph = rd->redistribute(*initgraph_manip, true);

  return rcp_balanced_graph;
}


/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> FSI::BlockMonolithic::SwitchDomains(Teuchos::RCP<Epetra_Map> rownodes,
    std::map<int, int>* nodeOwner, Teuchos::RCP<Epetra_CrsGraph> bal_graph, const Epetra_Comm& comm)
{
  /******************/
  /* switch domains */
  /******************/
  /*
   * The distribution of the nodes to procs is done. Now we want to
   * switch these domains of nodes among the procs (if necessary) such
   * that at the interface fluid and structure nodes are handled by the
   * same proc.
   */

  interfaceprocs_.remove(-1);

  int numproc = comm.NumProc();
  int myrank = comm.MyPID();

  int maxNumOfNodes;
  int nN = rownodes->NumMyElements();
  comm.MaxAll(&nN, &maxNumOfNodes, 1);

  int numMyFutureNodes = rownodes->NumMyElements();
  int myFutureNodeGIDs[maxNumOfNodes];  // contains GIDs of proc's nodes, later needed for map
                                        // construction
  int* nodeIdsPtr = rownodes->MyGlobalElements();

  for (int k = 0; k < numMyFutureNodes;
       ++k)  // Just insert the old node IDs. Maybe not every proc will change something here.
    myFutureNodeGIDs[k] = nodeIdsPtr[k];

  int countNodes[numproc];  // We use this array to determine to which proc another proc has to send
                            // its nodes.

  double nonmatch = 0;  // How many of my nodes don't want to be with me?
  double ratio;         // nonmatch is double because later this ratio is needed.

  double numInterfaceNodes = 0;
  int globalInterfaceNodes = 0;

  int interfaceproc = 0;  // Does this proc own interface nodes?

  bool checkProcAgain;  // Always the current nodes of the proc are checked. If proc 0
                        // switches with proc 1, we have to check proc 0 again (such that
                        // the old nodes of proc 1 are checked). If we continue to proc 1,
                        // the old nodes of proc 0 would be checked again.

  int sendToProc = -1;
  int numInterfaceNodesCopy = 0;


  for (int proc = 0; proc < numproc; ++proc)
  {
    checkProcAgain = true;
    while (checkProcAgain)
    {
      checkProcAgain = false;
      numInterfaceNodes = 0;
      nonmatch = 0;
      ratio = 0.0;
      for (int i = 0; i < numproc; ++i)
      {
        countNodes[i] = 0;
      }

      int numNodes = numMyFutureNodes;
      int increase;
      int increaseCopy;
      int wishProc;  // the proc a node wants to belong to
      int wishProcCopy;
      int nodeIndex;

      comm.Broadcast(&numNodes, 1, proc);
      for (int i = 0; i < numNodes; ++i)
      {
        if (myrank == proc)
        {
          nodeIndex = myFutureNodeGIDs[i];
        }
        comm.Broadcast(&nodeIndex, 1, proc);
        increase = -1;
        wishProc = -1;
        try
        {                                       // Check if node is interface node. If yes:
          wishProc = nodeOwner->at(nodeIndex);  // Save procID to which node wants to belong.
          if (wishProc != proc)
          {
            increase = 1;
          }  // Node is NOT on the proc where it needs to be.
          else
          {
            increase = 0;
          }  // Node is on the proc where it belongs.
        }
        catch (std::exception& exc)
        {
        }
        for (int sender = 0; sender < numproc; ++sender)
        {
          increaseCopy = increase;
          wishProcCopy = wishProc;
          comm.Broadcast(&increaseCopy, 1, sender);
          comm.Broadcast(&wishProcCopy, 1, sender);
          if (increaseCopy != -1)
          {
            if (myrank == proc)
            {
              nonmatch += increaseCopy;
              ++countNodes[wishProcCopy];  // Save the info to which proc this nodes wants to
                                           // belong. If the node has to be sent away it will be
                                           // sent to the proc with the highest number of
                                           // countNodes.
              ++numInterfaceNodes;
            }
          }
        }
      }  // for (int i=0; i<numNodes; ++i){

      // Check if I have to give away my nodes and if yes, to which proc.

      if (myrank == proc && numInterfaceNodes != 0)
      {
        ratio = nonmatch / numInterfaceNodes;
        int sendToProcCopy = sendToProc;
        if (ratio > 0.5)
        {  // I need to give my nodes to another processor.
          sendToProc = 0;
          for (int p = 1; p < numproc; ++p)
          {
            if (countNodes[p] > countNodes[sendToProc]) sendToProc = p;
          }
          if (sendToProc == sendToProcCopy)
          {  // Case: The proc wants to send back its nodes to
            if (numInterfaceNodes <= numInterfaceNodesCopy)  // the proc from which it received
                                                             // them. Give them to the proc with
              sendToProc = -1;  // the biggest share of the opposing patch. Explanation why we need
                                // this below.
          }
          else
            interfaceproc = 0;
        }
        else
        {
          sendToProc = -1;  // I can keep my nodes.
          interfaceproc = 1;
        }
      }
      else
      {
        sendToProc = -1;
      }

      /*
       * Explanation of the case: Proc wants to send back its nodes to the proc from which it
       * received them
       *
       * Some problems are sensible to the selected edge weights in the graph which results in a bad
       * partition where maybe one large patch on the undistributed domain is opposed to two small
       * patches which together match the big domain. Without detecting this case this results in an
       * infinite loop sending the nodes from one proc to the other and back. Current
       * implementation: Check which proc has the larger share of the opposing patch. Still not a
       * good distribution but the best we can achieve with the given decomposition.
       */

      numInterfaceNodesCopy =
          numInterfaceNodes;  // Save this value in the copy-variable until next loop iteration.

      //        if (comm_.MyPID()==proc)
      //          std::cout<<"\nProc "<<proc<<" has to switch with "<<sendToProc<<", it is ratio =
      //          "<<ratio<<std::endl;

      comm.Broadcast(&sendToProc, 1, proc);

      // Now do the actual switch between two procs, proc and sendToProc.

      int numNodesForwards;   // Number of nodes given from proc to sendToProc.
      int numNodesBackwards;  // Number of nodes given from sendToProc to proc.

      int nodeForwards;   // Nodes given from proc to sendToProc.
      int nodeBackwards;  // Nodes given from sendToProc to proc.

      if (sendToProc != -1 && sendToProc != proc)
      {
        for (int i = 0; i < maxNumOfNodes; ++i)
        {
          if (myrank == proc)
          {
            nodeForwards = myFutureNodeGIDs[i];
          }
          comm.Broadcast(&nodeForwards, 1, proc);
          if (myrank == sendToProc)
          {
            nodeBackwards = myFutureNodeGIDs[i];
            myFutureNodeGIDs[i] = nodeForwards;
          }
          comm.Broadcast(&nodeBackwards, 1, sendToProc);
          if (myrank == proc) myFutureNodeGIDs[i] = nodeBackwards;
        }
        if (myrank == proc)
        {
          numNodesForwards = numMyFutureNodes;
        }
        comm.Broadcast(&numNodesForwards, 1, proc);
        if (myrank == sendToProc)
        {
          numNodesBackwards = numMyFutureNodes;
          numMyFutureNodes = numNodesForwards;
          interfaceproc = 1;
        }
        comm.Broadcast(&numNodesBackwards, 1, sendToProc);
        if (myrank == proc) numMyFutureNodes = numNodesBackwards;

        checkProcAgain = true;  // Procs switched nodes. Now let the same proc check its new nodes.
      }

      comm.Broadcast(&numInterfaceNodes, 1, proc);
      globalInterfaceNodes += numInterfaceNodes;
    }
  }  // for (int proc=0; proc<numproc; ++proc){

  // Which procs own interface nodes?
  for (int p = 0; p < numproc; ++p)
  {
    int interfaceproccopy = interfaceproc;
    comm.Broadcast(&interfaceproccopy, 1, p);
    if (interfaceproccopy) interfaceprocs_.push_back(p);
  }

  //    if (myrank == 0){
  //      std::list<int>::iterator listit;
  //      for (listit = interfaceprocs_.begin(); listit != interfaceprocs_.end(); ++listit)
  //        std::cout<<"\nInterfaceproc: "<<*listit;
  //    }
  //    comm_.Barrier();

  // Create new map, new graph, an Epetra Exporter and do the new distribution.
  Epetra_Map newDistribution(-1, numMyFutureNodes, myFutureNodeGIDs, 0, comm);

  Teuchos::RCP<Epetra_CrsGraph> switched_bal_graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, newDistribution, 0));
  Epetra_Export exporter(*rownodes, newDistribution);

  int switch_err = switched_bal_graph->Export(*bal_graph, exporter, Insert);

  if (switch_err) dserror("Switch of domains not successful, error code %d!", switch_err);

  return switched_bal_graph;
}

/*---------------------------------------------------------------------------*/
void FSI::BlockMonolithic::RestoreRedistStructFluidGraph(
    std::map<int, std::vector<int>>& edgesToRemove, std::map<int, std::vector<int>>& edgesToInsert,
    Teuchos::RCP<Epetra_CrsGraph> monolithicGraph, Teuchos::RCP<Epetra_Map> monolithicRowmap,
    Teuchos::RCP<Epetra_Map> monolithicColmap, Teuchos::RCP<Epetra_CrsGraph> structureGraphRedist,
    Teuchos::RCP<Epetra_CrsGraph> fluidGraphRedist,
    std::map<int, std::vector<int>>& fluidToStructureMap)
{
  // todo: edgesToRemove not necessarily needed. just remove every fluid col node from structure row
  // and vice versa with info from what is a fluid node

  Epetra_Map finalStructureMap = (Epetra_Map&)structureGraphRedist->RowMap();
  Epetra_Map finalFluidMap = (Epetra_Map&)fluidGraphRedist->RowMap();

  int numMyNodes = monolithicRowmap->NumMyElements();
  int maxNumIndices = monolithicGraph->MaxNumIndices();
  int actNumMyNodes = 0;
  int ind[maxNumIndices];

  for (int i = 0; i < numMyNodes; ++i)
  {
    int gid = monolithicRowmap->GID(i);
    int err = monolithicGraph->ExtractGlobalRowCopy(gid, maxNumIndices, actNumMyNodes, ind);
    if (err != 0) dserror("\nExtractGlobalRowCopy failed, error code %d!", err);

    std::vector<int> insertCoupl;
    std::vector<int> removeCoupl;

    // edges to insert?
    int numInsert = 0;
    try
    {
      // yes
      insertCoupl = edgesToInsert.at(gid);
      numInsert = insertCoupl.size();
    }
    catch (std::exception& exc)
    {
    }  // no

    // edges to remove?
    int numRemove = 0;
    try
    {
      // yes
      removeCoupl = edgesToRemove.at(gid);
      numRemove = removeCoupl.size();
    }
    catch (std::exception& exc)
    {
    }  // no

    std::vector<int> structureInd;
    std::vector<int> fluidInd;

    // First deal with nodes that were already present in monolithicGraph. Just don't insert nodes
    // included in removeCoupl

    int countStructure = 0;
    int countFluid = 0;

    // check if fluid or structure ROW

    int rowlid = finalStructureMap.LID(gid);
    if (rowlid == -1)
    {  // not a structure row node
      rowlid = finalFluidMap.LID(gid);
      if (rowlid == -1)  // not a fluid row node either
        dserror("\nNode has to be either structure or fluid node!");
      else  // fluid row node
      {
        for (int j = 0; j < actNumMyNodes; ++j)
        {
          if (numRemove > 0)
          {
            bool insert = true;
            for (int k = 0; k < numRemove; ++k)
            {
              if (removeCoupl[k] == ind[j])
              {
                insert = false;
                break;
              }
            }
            if (insert == true)
            {
              // edge will be inserted
              fluidInd.push_back(ind[j]);
              countFluid++;
            }
          }  // if (numRemove>0)

          else  // nothing to remove
          {
            fluidInd.push_back(ind[j]);
            countFluid++;
          }
        }

        // insertion

        if (numInsert > 0)
        {
          for (int k = 0; k < numInsert; ++k)
          {
            fluidInd.push_back(insertCoupl[k]);
            countFluid++;
          }
        }
        // else: do nothing
      }
    }

    else  // structure row node
    {
      for (int j = 0; j < actNumMyNodes; ++j)
      {
        if (numRemove > 0)
        {
          bool insert = true;
          for (int k = 0; k < numRemove; ++k)
          {
            if (removeCoupl[k] == ind[j])
            {
              insert = false;
              break;
            }
          }
          if (insert == true)
          {
            // edge will be inserted
            structureInd.push_back(ind[j]);
            countStructure++;
          }
        }     // if (numRemove>0)
        else  // nothing to remove
        {
          structureInd.push_back(ind[j]);
          countStructure++;
        }
      }

      // insertion

      // Not needed here. No structure nodes have been removed previously.
    }


    if (countFluid > 0)
    {
      int err = fluidGraphRedist->InsertGlobalIndices(gid, countFluid, &fluidInd[0]);
      if (err != 0) dserror("\nInsert global indices failed, error code %d!", err);
    }
    else if (countStructure > 0)
    {
      int err = structureGraphRedist->InsertGlobalIndices(gid, countStructure, &structureInd[0]);
      if (err != 0) dserror("\nInsert global indices failed, error code %d!", err);
    }
  }


  // FINALLY: insert the deleted rows connected to fluid interface nodes

  std::map<int, std::vector<int>>::iterator it;

  for (it = fluidToStructureMap.begin(); it != fluidToStructureMap.end(); it++)
  {
    int gid = it->first;
    int lid = finalFluidMap.LID(gid);

    if (lid != -1)
    {
      std::vector<int> insert = edgesToInsert[gid];
      int num = insert.size();

      int err = fluidGraphRedist->InsertGlobalIndices(gid, num, &insert[0]);
      if (err != 0) dserror("\nInsert global indices failed, error code %d!", err);
    }
  }
}

/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::InsertDeletedEdges(std::map<int, std::list<int>>* deletedEdges,
    Teuchos::RCP<Epetra_Map> switched_rownodes, Teuchos::RCP<Epetra_CrsGraph> switched_bal_graph)
{
  // Insert deleted edges
  std::map<int, std::list<int>>::iterator delEdgeIt;
  int row;
  std::list<int> col_list;
  std::list<int>::iterator col_list_it;
  for (delEdgeIt = deletedEdges->begin(); delEdgeIt != deletedEdges->end(); delEdgeIt++)
  {
    row = delEdgeIt->first;
    if (switched_rownodes->LID(row) != -1)
    {
      col_list = delEdgeIt->second;
      int col_length = col_list.size();
      int col_indices[col_length];
      int k = 0;
      for (col_list_it = col_list.begin(); col_list_it != col_list.end(); col_list_it++)
      {
        col_indices[k] = *col_list_it;
        ++k;
      }
      int success = switched_bal_graph->InsertGlobalIndices(row, col_length, col_indices);
      if (success != 0) dserror("\nInsert global indices failed, error code %d!", success);
    }
  }

  switched_bal_graph->FillComplete();
  switched_bal_graph->OptimizeStorage();
}

/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::FindNodeRelatedToDof(std::map<int, DRT::Node*>* nodes, int gdofid,
    Teuchos::RCP<DRT::Discretization> discretization, int* re)
{
  re[0] = -2;  // code: the node cannot be found on this proc
  bool breakout = false;
  std::vector<int> dofs;
  std::map<int, DRT::Node*>::iterator nodeiterator;

  for (nodeiterator = nodes->begin(); nodeiterator != nodes->end(); nodeiterator++)
  {
    discretization->Dof((const DRT::Node*)nodeiterator->second, dofs);
    for (int i = 0; i < (int)dofs.size(); ++i)
    {
      if (dofs[i] == gdofid)
      {
        re[0] = nodeiterator->first;
        re[1] = nodeiterator->second->Owner();
        breakout = true;
        break;
      }
    }
    if (breakout == true) break;
    dofs.clear();
  }
}

void FSI::BlockMonolithic::BuildMonolithicGraph(Teuchos::RCP<Epetra_CrsGraph> monolithicGraph,
    std::map<int, std::vector<int>>& deletedEdges, std::map<int, std::vector<int>>& insertedEdges,
    std::map<int, std::vector<int>>& fluidToStructureMap,
    std::map<int, std::vector<int>>& structureToFluidMap,
    Teuchos::RCP<DRT::Discretization> structuredis,  ///< structure discretization
    Teuchos::RCP<DRT::Discretization> fluiddis)
{
  int numproc = monolithicGraph->Comm().NumProc();
  int myrank = monolithicGraph->Comm().MyPID();

  // create nodal graphs and maps of structure and fluid field

  Teuchos::RCP<const Epetra_CrsGraph> structureGraph = structuredis->BuildNodeGraph();
  Teuchos::RCP<const Epetra_CrsGraph> fluidGraph = fluiddis->BuildNodeGraph();
  const Epetra_Map& structureGraph_map = (Epetra_Map&)structureGraph->RowMap();
  const Epetra_Map& fluidGraph_map = (Epetra_Map&)fluidGraph->RowMap();
  //
  //  // hypergraph edge weights: set all values to weight1
  //  crs_hge_weights->PutScalar(weight1);
  //  double weight = weight2;
  //

  /*******************************************/
  /* copy line after line in monolithicGraph */
  /*******************************************/
  // Note: Fluid interface lines are omitted. Respective Columns are omitted too.


  // fluid part: nothing is inserted in monolithicGraph, only information is collected
  // insertion happens later

  int numMyFluidNodes = fluidGraph_map.NumMyElements();
  int maxNumFluidIndices = fluidGraph->MaxNumIndices();
  int actNumIndices = 0;  // actual number of entries in "indices"
  int indices_f[maxNumFluidIndices];
  std::vector<int> coupling;
  bool insertCoupling = false;

  std::map<int, std::vector<int>>
      strNodesCouplToFluid;  // which structure interface node is coupled to which fluid nodes
  std::map<int, std::vector<int>>
      fluidInterfNodesConnect;  // which fluid interface node is coupled to which fluid nodes
  std::map<int, std::vector<int>> fluidFieldNodesConnect;  // connectivity of fluid field nodes

  for (int i = 0; i < numMyFluidNodes; ++i)
  {
    int gid = fluidGraph_map.GID(i);

    // extract indices

    int err = fluidGraph->ExtractGlobalRowCopy(gid, maxNumFluidIndices, actNumIndices, indices_f);
    if (err != 0) dserror("ExtractGlobalRowCopy failed, error = %d!", err);

    // check if interface node

    try
    {
      coupling = fluidToStructureMap.at(gid);  // Try to get vector with coupling entries
      insertCoupling = true;
      // crs_hge_weights->ReplaceGlobalValues(1,&weight,&gid);
    }
    catch (std::exception& exc)
    {
    }  // Not an interface node. Do nothing


    // get insertion information

    if (insertCoupling == false)  // field node
    {
      // do not insert fluid interface nodes

      std::vector<int> insert_ind;
      int numInsert = 0;

      // save deleted edges
      for (int k = 0; k < actNumIndices; ++k)
      {
        try
        {
          fluidToStructureMap.at(indices_f[k]);
          deletedEdges[gid].push_back(indices_f[k]);
        }
        catch (std::exception& exc)
        {
          insert_ind.push_back(indices_f[k]);
          numInsert++;
        }
      }

      fluidFieldNodesConnect[gid] = insert_ind;
    }
    else  // interface node
    {
      // cast array to std::vector
      std::vector<int> indices_f_coupl(indices_f, indices_f + actNumIndices);

      for (int v = 0; v < actNumIndices; ++v)
      {
        try
        {
          fluidToStructureMap.at(indices_f_coupl[v]);  // don't consider fluid interface nodes
        }
        catch (std::exception& exc)
        {
          if (indices_f_coupl[v] != gid)
          {
            strNodesCouplToFluid[coupling[0]].push_back(
                indices_f_coupl[v]);  // structure nodes --- coupled fluid nodes
            fluidInterfNodesConnect[gid].push_back(
                indices_f_coupl[v]);  // fluid nodes ----- coupled fluid nodes
          }
          fluidInterfNodesConnect[gid].push_back(gid);
        }
      }

      // save deleted edges
      for (int k = 0; k < actNumIndices; ++k)
      {
        deletedEdges[gid].push_back(indices_f[k]);
      }

      insertCoupling = false;
    }
  }  // end of fluid part I



  // communicate connectivity information of structure interface nodes among all procs

  std::map<int, std::vector<int>> connectInfo;  // structure interface node to fluid field nodes

  std::vector<int> transferVec;
  int transfer[1000];
  std::map<int, std::vector<int>>::iterator mapIter;

  for (int p = 0; p < numproc; ++p)
  {
    int numNodes = strNodesCouplToFluid.size();
    monolithicGraph->Comm().Broadcast(&numNodes, 1, p);

    mapIter = strNodesCouplToFluid.begin();

    for (int n = 0; n < numNodes; ++n)
    {
      int gid = 0;
      int numCoupl = 0;

      if (myrank == p)
      {
        gid = mapIter->first;
        numCoupl = (int)strNodesCouplToFluid[gid].size();
        transferVec = mapIter->second;
        for (int t = 0; t < numCoupl; ++t) transfer[t] = transferVec[t];
      }

      monolithicGraph->Comm().Broadcast(&gid, 1, p);
      monolithicGraph->Comm().Broadcast(&numCoupl, 1, p);
      monolithicGraph->Comm().Broadcast(transfer, numCoupl, p);

      std::vector<int> nodes;

      for (int a = 0; a < numCoupl; ++a)
      {
        nodes.push_back(transfer[a]);
      }

      connectInfo[gid] = nodes;

      if (myrank == p) mapIter++;
    }

    // communicate deleted edges

    numNodes = deletedEdges.size();
    monolithicGraph->Comm().Broadcast(&numNodes, 1, p);

    mapIter = deletedEdges.begin();

    for (int n = 0; n < numNodes; ++n)
    {
      int gid = 0;
      int numDelEd = 0;

      if (myrank == p)
      {
        gid = mapIter->first;
        numDelEd = (int)deletedEdges[gid].size();
        transferVec = mapIter->second;
        for (int t = 0; t < numDelEd; ++t) transfer[t] = transferVec[t];
      }

      monolithicGraph->Comm().Broadcast(&gid, 1, p);
      monolithicGraph->Comm().Broadcast(&numDelEd, 1, p);
      monolithicGraph->Comm().Broadcast(transfer, numDelEd, p);

      std::vector<int> delEd;

      for (int a = 0; a < numDelEd; ++a)
      {
        delEd.push_back(transfer[a]);
      }

      deletedEdges[gid] = delEd;

      if (myrank == p) mapIter++;
    }


  }  // end of communication


  // build inverse map

  std::map<int, std::vector<int>>
      inverseConnectInfo;  // fluid field nodes to structure interface nodes
  std::map<int, std::vector<int>>::iterator connectIt;

  for (connectIt = connectInfo.begin(); connectIt != connectInfo.end(); connectIt++)
  {
    int first = connectIt->first;
    std::vector<int> second = connectIt->second;
    int num = second.size();

    for (int i = 0; i < num; ++i)
    {
      try
      {
        inverseConnectInfo[second[i]].at(
            first);  // check if "first" is already contained in inverse map
      }
      catch (std::exception& exc)
      {
        inverseConnectInfo[second[i]].push_back(first);
      }
    }
  }

  // Merge both maps. This is needed later for removal of inserted edges.

  insertedEdges = connectInfo;
  insertedEdges.insert(inverseConnectInfo.begin(), inverseConnectInfo.end());

  // end of map building

  // structure part

  int numMyStructureNodes = structureGraph_map.NumMyElements();
  int maxNumStructureIndices = structureGraph->MaxNumIndices();
  int indices_s[maxNumStructureIndices];

  for (int i = 0; i < numMyStructureNodes; ++i)
  {
    int gid = structureGraph_map.GID(i);

    // extract indices

    int err =
        structureGraph->ExtractGlobalRowCopy(gid, maxNumStructureIndices, actNumIndices, indices_s);
    if (err != 0) dserror("ExtractGlobalRowCopy failed, error = %d!", err);


    // Check if interface node. If so: insert coupling entries from fluid

    try
    {
      coupling = structureToFluidMap.at(gid);  // Try to get vector with coupling entries
      insertCoupling = true;
      // crs_hge_weights->ReplaceGlobalValues(1,&weight,&gid);
    }
    catch (std::exception& exc)
    {
    }  // Not an interface node. Do nothing


    if (insertCoupling == false)
    {
      // do not insert indices at fluid interface

      std::vector<int> insert_ind;
      int numInsert = 0;

      // save deleted edges
      for (int k = 0; k < actNumIndices; ++k)
      {
        try
        {
          fluidToStructureMap.at(indices_s[k]);
          deletedEdges[gid].push_back(indices_s[k]);
        }
        catch (std::exception& exc)
        {
          insert_ind.push_back(indices_s[k]);
          numInsert++;
        }
      }

      err = monolithicGraph->InsertGlobalIndices(gid, numInsert, &insert_ind[0]);
      if (err != 0) dserror("InsertGlobalIndices failed, error = %d!", err);
    }
    else
    {
      // do not insert indices at fluid interface

      std::vector<int> insert_ind;
      int numInsert = 0;
      std::vector<int> delEd = deletedEdges[gid];

      // save deleted edges
      for (int k = 0; k < actNumIndices; ++k)
      {
        try
        {
          fluidToStructureMap.at(indices_s[k]);
          delEd.push_back(indices_s[k]);
        }
        catch (std::exception& exc)
        {
          insert_ind.push_back(indices_s[k]);
          numInsert++;
        }
      }

      std::vector<int> nodesFromFluid = connectInfo[gid];

      int numTransfer = nodesFromFluid.size();
      for (int k = 0; k < numTransfer; ++k)
      {
        try
        {
          fluidToStructureMap.at(nodesFromFluid[k]);
          deletedEdges[gid].push_back(nodesFromFluid[k]);
        }
        catch (std::exception& exc)
        {
          insert_ind.push_back(nodesFromFluid[k]);
          numInsert++;
        }
      }

      numInsert = insert_ind.size();

      err = monolithicGraph->InsertGlobalIndices(gid, numInsert, &insert_ind[0]);
      if (err != 0) dserror("InsertGlobalIndices failed, error = %d!", err);
    }
  }  // end of structure part



  // now run again through fluid rows to insert elements for symmetric graph

  // for (connectIt = inverseConnectInfo.begin(); connectIt != inverseConnectInfo.end();
  // ++connectIt)
  for (int i = 0; i < numMyFluidNodes; ++i)
  {
    int gid = fluidGraph_map.GID(i);
    // int lid = fluidGraph_map.LID(gid);

    //    if (lid != -1)
    //    {

    try
    {
      fluidToStructureMap.at(gid);  // Interface node. Do nothing
    }
    catch (std::exception& exc)
    {
      std::vector<int> insert_ind;
      int numInsert = 0;

      // coupling nodes

      try
      {
        std::vector<int> nodesFromStructure = inverseConnectInfo.at(gid);

        // std::vector<int> nodesFromStructure = inverseConnectInfo[gid];

        int numTransfer = nodesFromStructure.size();
        for (int k = 0; k < numTransfer; ++k)
        {
          try
          {
            fluidToStructureMap.at(nodesFromStructure[k]);
            deletedEdges[gid].push_back(nodesFromStructure[k]);
          }
          catch (std::exception& exc)
          {
            insert_ind.push_back(nodesFromStructure[k]);
            numInsert++;
          }
        }
      }
      catch (std::exception& exc)
      {
      }

      // field nodes

      insert_ind.insert(
          insert_ind.end(), fluidFieldNodesConnect[gid].begin(), fluidFieldNodesConnect[gid].end());
      numInsert += (int)fluidFieldNodesConnect[gid].size();

      //        std::cout<<"gid: "<<gid<<" rank: "<<myrank<<" numInsert: "<<numInsert<<std::endl;
      //        for (int m=0; m<numInsert; ++m)
      //          std::cout<<"rank: "<<myrank<<" ind: "<<insert_ind[m]<<std::endl;

      // insertion

      int err = monolithicGraph->InsertGlobalIndices(gid, numInsert, &insert_ind[0]);
      if (err != 0) dserror("InsertGlobalIndices failed, error = %d!", err);
    }

    //}

  }  // end of fluid part II

  monolithicGraph->FillComplete();
}


// Teuchos::RCP<Epetra_Map> FSI::BlockMonolithic::GetStructureRowMap(
//    const Epetra_Map& oldStructureMap,
//    Teuchos::RCP<Epetra_Map> monolithicRownodes,
//    std::map<int,std::vector<int> > &fluidToStructureMap)
//{
//
//  int numproc = oldStructureMap.Comm().NumProc();
//  int myrank = oldStructureMap.Comm().MyPID();
//
//  // build arrays with all structure and all fluid nodes from all procs
//
//  int numGlobalStructureNodes = oldStructureMap.NumGlobalElements();
//  int numMyStructureNodes = oldStructureMap.NumMyElements();
//
//  std::vector<int> allStructureNodes;
//
//  int* myStructureNodes = oldStructureMap.MyGlobalElements();
//
//  // communicate information among procs
//
//  for (int p=0; p<numproc; ++p)
//  {
//    int numReceivedStructureNodes = numMyStructureNodes;
//    oldStructureMap.Comm().Broadcast(&numReceivedStructureNodes,1,p);
//
//    int receivedStructureNodes[numReceivedStructureNodes];
//
//    if (myrank == p)
//    {
//      for (int i=0; i<numReceivedStructureNodes; ++i)
//        receivedStructureNodes[i] = myStructureNodes[i];
//    }
//
//    oldStructureMap.Comm().Broadcast(receivedStructureNodes,numReceivedStructureNodes,p);
//
//    for (int i=0; i<numReceivedStructureNodes; ++i)
//      allStructureNodes.push_back(receivedStructureNodes[i]);
//
//  } // end of communication
//
//
//  // separate fluid and structure nodes from monolithic maps
//
//  std::vector<int> myRedistStructureRowNodes;
//  int numMyRedistStructureRowNodes = 0;
//
//
//  for (int i=0; i<numGlobalStructureNodes; ++i)
//  {
//    int lid = monolithicRownodes->LID(allStructureNodes[i]);
//    if (lid != -1)  // proc owns this structure node
//    {
//      myRedistStructureRowNodes.push_back(allStructureNodes[i]);
//      numMyRedistStructureRowNodes++;
//    }
//  }
//
//  Teuchos::RCP<Epetra_Map> structureRowmap = Teuchos::rcp(new
//  Epetra_Map(numGlobalStructureNodes,numMyRedistStructureRowNodes,&myRedistStructureRowNodes[0],0,oldStructureMap.Comm()));
//
//  return structureRowmap;
//
//}
//
//
// Teuchos::RCP<Epetra_Map> FSI::BlockMonolithic::GetFluidRowMap(
//    const Epetra_Map& oldFluidMap,
//    Teuchos::RCP<Epetra_Map> monolithicRownodes,
//    std::map<int,std::vector<int> > &fluidToStructureMap,
//    std::vector<int> &allFluidNodes)
//{
//
//  int numproc = oldFluidMap.Comm().NumProc();
//  int myrank = oldFluidMap.Comm().MyPID();
//
//  // build arrays with all structure and all fluid nodes from all procs
//
//  int numGlobalFluidNodes = oldFluidMap.NumGlobalElements();
//  int numMyFluidNodes = oldFluidMap.NumMyElements();
//
//  int* myFluidNodes = oldFluidMap.MyGlobalElements();
//
//  // communicate information among procs
//
//  for (int p=0; p<numproc; ++p)
//  {
//
//    int numReceivedFluidNodes = numMyFluidNodes;
//    oldFluidMap.Comm().Broadcast(&numReceivedFluidNodes,1,p);
//
//    int receivedFluidNodes[numReceivedFluidNodes];
//
//    if (myrank == p)
//    {
//      for (int i=0; i<numReceivedFluidNodes; ++i)
//        receivedFluidNodes[i] = myFluidNodes[i];
//    }
//
//    oldFluidMap.Comm().Broadcast(receivedFluidNodes,numReceivedFluidNodes,p);
//
//    for (int i=0; i<numReceivedFluidNodes; ++i)
//      allFluidNodes.push_back(receivedFluidNodes[i]);
//
//  } // end of communication
//
//
//  // separate fluid and structure nodes from monolithic maps
//
//  std::vector<int> myRedistFluidRowNodes;
//  int numMyRedistFluidRowNodes = 0;
//
//  for (int i=0; i<numGlobalFluidNodes; ++i)
//  {
//    int lid = monolithicRownodes->LID(allFluidNodes[i]);
//    if (lid != -1)  // proc owns this structure node
//    {
//      try
//      {
//        fluidToStructureMap.at(allFluidNodes[i]);
//      }
//      catch (std::exception& exc)
//      {
//        myRedistFluidRowNodes.push_back(allFluidNodes[i]);
//        numMyRedistFluidRowNodes++;
//      }
//    }
//  }
//
//
//  // fluid: in monolithicGraph all fluid interface nodes have been omitted. Insert now
//
//  std::map<int, std::vector<int> >::iterator it;
//
//  for (it = fluidToStructureMap.begin(); it != fluidToStructureMap.end(); ++it)
//  {
//
//    int structGID = it->second[0];
//    int structLID = monolithicRownodes->LID(structGID);
//
//    if (structLID != -1)
//    {
//      int fluidGID = it->first;
//      myRedistFluidRowNodes.push_back(fluidGID);
//      numMyRedistFluidRowNodes++;
//    }
//  }
//
//  Teuchos::RCP<Epetra_Map> fluidRowmap = Teuchos::rcp(new
//  Epetra_Map(numGlobalFluidNodes,numMyRedistFluidRowNodes,&myRedistFluidRowNodes[0],0,oldFluidMap.Comm()));
//
//  return fluidRowmap;
//
//}


Teuchos::RCP<Epetra_Map> FSI::BlockMonolithic::GetRedistRowMap(const Epetra_Map& oldMap,
    Teuchos::RCP<Epetra_Map> monolithicRownodes,
    std::map<int, std::vector<int>>& fluidToStructureMap, bool fluid)
{
  int numproc = oldMap.Comm().NumProc();
  int myrank = oldMap.Comm().MyPID();

  // build arrays with all structure and all fluid nodes from all procs

  int numGlobalNodes = oldMap.NumGlobalElements();
  int numMyNodes = oldMap.NumMyElements();

  int* myNodes = oldMap.MyGlobalElements();

  std::vector<int> allNodes;

  // communicate information among procs

  for (int p = 0; p < numproc; ++p)
  {
    int numReceivedNodes = numMyNodes;
    oldMap.Comm().Broadcast(&numReceivedNodes, 1, p);

    int receivedNodes[numReceivedNodes];

    if (myrank == p)
    {
      for (int i = 0; i < numReceivedNodes; ++i) receivedNodes[i] = myNodes[i];
    }

    oldMap.Comm().Broadcast(receivedNodes, numReceivedNodes, p);

    for (int i = 0; i < numReceivedNodes; ++i) allNodes.push_back(receivedNodes[i]);

  }  // end of communication


  // separate fluid and structure nodes from monolithic maps

  std::vector<int> myRedistRowNodes;
  int numMyRedistRowNodes = 0;

  for (int i = 0; i < numGlobalNodes; ++i)
  {
    int lid = monolithicRownodes->LID(allNodes[i]);
    if (lid != -1)  // proc owns this structure node
    {
      try
      {
        fluidToStructureMap.at(allNodes[i]);
      }
      catch (std::exception& exc)
      {
        myRedistRowNodes.push_back(allNodes[i]);
        numMyRedistRowNodes++;
      }
    }
  }


  // fluid: in monolithicGraph all fluid interface nodes have been omitted. Insert now
  if (fluid == true)
  {
    std::map<int, std::vector<int>>::iterator it;

    for (it = fluidToStructureMap.begin(); it != fluidToStructureMap.end(); ++it)
    {
      int structGID = it->second[0];
      int structLID = monolithicRownodes->LID(structGID);

      if (structLID != -1)
      {
        int fluidGID = it->first;
        myRedistRowNodes.push_back(fluidGID);
        numMyRedistRowNodes++;
      }
    }
  }

  Teuchos::RCP<Epetra_Map> redistRowmap = Teuchos::rcp(
      new Epetra_Map(numGlobalNodes, numMyRedistRowNodes, &myRedistRowNodes[0], 0, oldMap.Comm()));

  return redistRowmap;
}
