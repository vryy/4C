/*----------------------------------------------------------------------------*/
/*!
\file fsi_monolithic_redistribution.cpp

\brief Parallel Domain Redistribution of monolithic FSI

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
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
void FSI::BlockMonolithic::RedistributeDomainDecomposition(
    const INPAR::FSI::Redistribute domain,
    const FSI_COUPLING coupling, const double inputWeight1,
    const double inputWeight2, const Epetra_Comm& comm)
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

  //initialize maps for elements
  std::map<int, Teuchos::RCP<DRT::Element> > structureelements;
  std::map<int, Teuchos::RCP<DRT::Element> > fluidelements;


  // access the discretizations
  Teuchos::RCP<DRT::Discretization> structuredis = DRT::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis = DRT::Problem::Instance()->GetDis("ale");

  // Fill maps based on condition for master side (masterdis != slavedis)
  DRT::UTILS::FindConditionObjects(*structuredis, structurenodes, structuregnodes,
      structureelements, "FSICoupling");

  // Fill maps based on condition for slave side (masterdis != slavedis)
  DRT::UTILS::FindConditionObjects(*fluiddis, fluidnodes, fluidgnodes,
      fluidelements, "FSICoupling");

  std::map<int, DRT::Node*>* slavenodesPtr = NULL;
  std::map<int, DRT::Node*>* mastergnodesPtr = NULL;

  if (coupling == fsi_iter_mortar_monolithicfluidsplit or coupling == fsi_iter_sliding_monolithicfluidsplit){
    slavenodesPtr = &fluidnodes;
    mastergnodesPtr = &structuregnodes;
  }
  else if (coupling == fsi_iter_mortar_monolithicstructuresplit or coupling == fsi_iter_sliding_monolithicstructuresplit){
    slavenodesPtr = &structurenodes;
    mastergnodesPtr = &fluidgnodes;
  }
  else
    dserror("\nDomain redistribution / rebalancing only implemented for mortar coupling!");


  /*******************************************/
  /* distribute masternodes to future owners */
  /*******************************************/

  std::map<int,int> nodeOwner;    // maps global node id and the owner of the neighbouring node
  std::map<int,int>* nodeOwnerPtr = &nodeOwner;
  std::map<int,std::list<int> > inverseNodeOwner;   // Create a map that maps Owner <-> Nodes
  std::map<int,std::list<int> >* inverseNodeOwnerPtr = &inverseNodeOwner;

  CreateNodeOwnerRelationship(nodeOwnerPtr, inverseNodeOwnerPtr, slavenodesPtr,
      mastergnodesPtr, structuredis, fluiddis, domain);

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
  const Epetra_Map& initgraph_map = (Epetra_Map&) initgraph->RowMap();

  Teuchos::RCP<Epetra_CrsMatrix> crs_ge_weights = Teuchos::rcp(new Epetra_CrsMatrix(Copy, initgraph_map, 0));
  Teuchos::RCP<Epetra_CrsGraph> initgraph_manip = Teuchos::rcp(new Epetra_CrsGraph(Copy, initgraph_map, 0));


  std::map<int,std::list<int> > deletedEdges;
  std::map<int,std::list<int> >* deletedEdgesPtr = &deletedEdges;

  BuildWeightedGraph(crs_ge_weights, initgraph_manip, initgraph, inputWeight1,
      inputWeight2, nodeOwnerPtr, inverseNodeOwnerPtr, deletedEdgesPtr, comm);

  /******************/
  /* redistribution */
  /******************/

  Teuchos::RCP<Epetra_CrsGraph> bal_graph = CallZoltan(initgraph_manip, crs_ge_weights);

  // Extract row map
  Teuchos::RCP<Epetra_Map> rownodes;
  rownodes = Teuchos::rcp(
      new Epetra_Map(-1, bal_graph->RowMap().NumMyElements(),
          bal_graph->RowMap().MyGlobalElements(), 0, comm));

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

  // Extract row map from the final graph. Column map will be done in a while, after we insert the deleted edges from above.
  Teuchos::RCP<Epetra_Map> switched_rownodes;
  switched_rownodes = Teuchos::rcp(
        new Epetra_Map(-1, switched_bal_graph->RowMap().NumMyElements(),
            switched_bal_graph->RowMap().MyGlobalElements(), 0, comm));


  // TODO: get only rowmap from switchdomains and don't build graph.

  /*
   * We inserted edges between all interface nodes on one proc in BuildWeightedGraph
   * which would lead to an incorrect column map if we used switched_bal_graph
   * to extract the column map. Therefore, we export the initgraph to the final
   * switched_rownodes and get the columnmap afterwards.
   */


  Teuchos::RCP<Epetra_CrsGraph> final_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*switched_rownodes,0));
  Epetra_Export exporter(initgraph_map, *switched_rownodes);

  int err = final_graph->Export(*initgraph,exporter,Insert);
    if (err)
      dserror("Export not successful, error code %d!",err);

  InsertDeletedEdges(deletedEdgesPtr, switched_rownodes, final_graph);

  // Now extract column map from the final graph.
  Teuchos::RCP<Epetra_Map> switched_colnodes;
  switched_colnodes = Teuchos::rcp(
      new Epetra_Map(-1, final_graph->ColMap().NumMyElements(),
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

  if (domain == INPAR::FSI::Redistribute_fluid){
    aledis->Redistribute(*switched_rownodes, *switched_colnodes);
  }

  // Distribute dofs in time integration vectors to right procs by creating time integrators new.
  // The Control File has to be rewritten as well.
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Teuchos::ParameterList> ioflags
    = Teuchos::rcp(new Teuchos::ParameterList(problem->IOParams()));
  Teuchos::RCP<IO::DiscretizationWriter> structureoutput = structuredis->Writer();
  Teuchos::RCP<IO::DiscretizationWriter> fluidoutput = fluiddis->Writer();
  Teuchos::RCP<IO::DiscretizationWriter> aleoutput = aledis->Writer();

  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();

  if (domain == INPAR::FSI::Redistribute_structure){
    fluidoutput->OverwriteResultFile();
    aleoutput->OverwriteResultFile();
    structureoutput->OverwriteResultFile();
    CreateStructureTimeIntegrator(fsidyn, structuredis);
    SetLambda();
  }
  else if (domain == INPAR::FSI::Redistribute_fluid){
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

  if (myrank==0)
    printf("Redistribution of domain in %f seconds.\n",timer.ElapsedTime());

  comm.Barrier();

  return;
}

/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::BuildWeightedGraph(
    Teuchos::RCP<Epetra_CrsMatrix> crs_ge_weights,
    Teuchos::RCP<Epetra_CrsGraph> initgraph_manip,
    Teuchos::RCP<const Epetra_CrsGraph> initgraph, const double inputWeight1,
    const double inputWeight2, std::map<int, int>* nodeOwner,
    std::map<int, std::list<int> >* inverseNodeOwner,
    std::map<int, std::list<int> >* deletedEdges, const Epetra_Comm& comm)
{
  /*********************************************/
  /* build nodal graph and insert edge weights */
  /*********************************************/

  int numproc = comm.NumProc();
  int myrank = comm.MyPID();

  const Epetra_Map& initgraph_map = (Epetra_Map&) initgraph->RowMap();

  int owner;        // owner in first loop
  int owner2;       // owner in second loop
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
  comm.MaxAll(&numNodes,&maxNodes,1);

  int numMyGraphRows = initgraph->NumMyRows();
  int maxNumGraphRows;
  comm.MaxAll(&numMyGraphRows,&maxNumGraphRows,1);

  for (int i=0; i<maxNumGraphRows; ++i){       // all procs have to be in the loop until the very last checked all its nodes

    int graphRowGID = crs_ge_weights->GRID(i);
    if (graphRowGID != -1){
      int success = initgraph->ExtractGlobalRowCopy(graphRowGID,graphLength,numGraphEntries,graphIndices);
      if (success != 0)
        dserror("\nExtract global values failed, error code %d!", success);
    }

    for (int proc=0; proc<numproc; ++proc){   // each processor checks if graphRowGID is an interface node

      int foundInterfaceNode = 0;
      noderowcopy = graphRowGID;
      comm.Broadcast(&noderowcopy,1,proc);
      numGraphEntriesCopy = numGraphEntries;
      comm.Broadcast(&numGraphEntriesCopy,1,proc);

      if (noderowcopy != -1){                 // only consider procs that still have graph rows in their map

        owner = -1;
        try{                                             // check which nodes are interface nodes
          owner = nodeOwner->at(noderowcopy);             // and if so: get (desired) owner
          foundInterfaceNode = 1;
          }
        catch (std::exception& exc)                      // not an interface node
        {
        }

        int foundInterfaceNodeSum;
        comm.SumAll(&foundInterfaceNode,&foundInterfaceNodeSum,1);

        if (myrank==proc){
          for (int c=0; c<numGraphEntries; ++c)
            graphIndicesCopy[c] = graphIndices[c];
        }

        comm.Broadcast(&graphIndicesCopy[0],numGraphEntriesCopy,proc);

        if (foundInterfaceNodeSum == 0){      // not an interface node
          if (initgraph_map.LID(noderowcopy) != -1){
            double fieldWeight[numGraphEntriesCopy];
            for (int w=0; w<numGraphEntriesCopy; ++w){
              fieldWeight[w] = inputWeight1;
            }

            int success;
            success = crs_ge_weights->InsertGlobalValues(noderowcopy,numGraphEntriesCopy,&fieldWeight[0],&graphIndicesCopy[0]);
            if (success != 0)
              dserror("\nInsert global values failed, error code %d!", success);
            success = initgraph_manip->InsertGlobalIndices(noderowcopy,numGraphEntriesCopy,&graphIndicesCopy[0]);
            if (success != 0)
              dserror("\nInsert global indices failed, error code %d!", success);

          }
        }
        else if(foundInterfaceNodeSum == 1){  // interface nodes

          for (int proc2=0; proc2<numproc; ++proc2){
            ownercopy = owner;
            comm.Broadcast(&ownercopy,1,proc2);

            if (ownercopy != -1){


              // we want to connect all interface nodes on one future proc with high weights
              std::list<int> interfaceNodesOnProc = (*inverseNodeOwner)[ownercopy];
              int listLength = interfaceNodesOnProc.size();

              // find interface node pairs

              double GIDweight;                                    // weight for graph position noderowcopy - graphIndices[k]
              double insertWeights[listLength];
              int insertIndices[listLength];

              int insert_k = 0;

              for (int k = 0; k < numGraphEntriesCopy; ++k){

                GIDweight = -1.0;
                int found = 0;
                owner2 = -1;

                if (graphIndicesCopy[k] != noderowcopy){

                  try{                                             // check which nodes are interface nodes
                    owner2 = nodeOwner->at(graphIndicesCopy[k]);   // and if so: get (desired) owner
                    if (owner2 == ownercopy){                      // node wants to be on the same proc as node related to noderowcopy
                      GIDweight = weight;
                    }
                    else{                                          // border between two neighboring patches
                      GIDweight = zero;
                    }
                    found = 1;
                  }
                  catch (std::exception& exc)                      // not an interface node
                  {
                  }

                }

                int foundSum;
                comm.SumAll(&found,&foundSum,1);

                if (foundSum == 1){                               // someone found an interface node

                  double GIDweightCopy;
                  for (int p=0; p<numproc; ++p){                  // collect weight info in insertWeights and insertIndices
                    GIDweightCopy = GIDweight;
                    comm.Broadcast(&GIDweightCopy,1,p);
                    if (GIDweightCopy != -1){
                      if (GIDweightCopy != zero){
                        insertWeights[insert_k] = GIDweightCopy;
                        insertIndices[insert_k] = graphIndicesCopy[k];
                        interfaceNodesOnProc.remove(graphIndicesCopy[k]);
                        insert_k++;
                      }
                      if (GIDweightCopy == zero){                 // delete that edge
                        //deletedEdges[noderowcopy] = graphIndicesCopy[k];
                        std::list<int>::iterator it = (*deletedEdges)[noderowcopy].begin();
                        (*deletedEdges)[noderowcopy].insert(it,graphIndicesCopy[k]);
                        interfaceNodesOnProc.remove(graphIndicesCopy[k]);
                      }
                      break;
                    }
                  }
                }
                else if (foundSum == 0){                           // no one found an interface node
                    insertWeights[insert_k] = inputWeight1;
                    insertIndices[insert_k] = graphIndicesCopy[k];
                    insert_k++;
                }
                else{                                              // we expect nodeOwner to be one to one, i.e. foundSum =  {0, 1}
                  dserror("\nNode <-> Owner mapping not one-to-one, %d procs have the same node in their maps!",foundSum);
                }

              }


              // connect interface nodes on one proc with high weights
              std::list<int>::iterator listIt = interfaceNodesOnProc.begin();
              listLength = interfaceNodesOnProc.size();

              for (int l=0; l<listLength; ++l){
                insertIndices[insert_k] = *listIt;
                insertWeights[insert_k] = weight;
                insert_k++;
                listIt++;
              }


              if (initgraph_map.LID(noderowcopy) != -1){           // get processor which accesses row noderowcopy and insert weights
                int success;
                success = crs_ge_weights->InsertGlobalValues(noderowcopy,insert_k,&insertWeights[0],&insertIndices[0]);
                if (success != 0)
                  dserror("\nReplace global values failed, error code %d!", success);
                success = initgraph_manip->InsertGlobalIndices(noderowcopy,insert_k,&insertIndices[0]);
                if (success != 0)
                  dserror("\nInsert global indices failed, error code %d!", success);
              }
            } // if (ownercopy != -1){
          } // for (int proc2=0; proc2<numproc; ++proc2){
        }
        else
          dserror("\n\nNode <-> Owner mapping not one-to-one, %d procs have the same node in their maps!",foundInterfaceNodeSum);

      } // if (noderowcopy != -1){

      comm.Barrier();

    } // for (int proc=0; proc<numproc; ++proc){

  } // for (int i=0; i<maxNodes; ++i){

  crs_ge_weights->FillComplete();
  crs_ge_weights->OptimizeStorage();

  initgraph_manip->FillComplete();
  initgraph_manip->OptimizeStorage();
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> FSI::BlockMonolithic::CallZoltan(
    Teuchos::RCP<const Epetra_CrsGraph> initgraph_manip,
    Teuchos::RCP<const Epetra_CrsMatrix> crs_ge_weights)
{

  Teuchos::ParameterList paramlist;

    paramlist.set("PARTITIONING METHOD", "GRAPH");
    paramlist.set("PRINT ZOLTAN METRICS", "2");
    Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
    sublist.set("GRAPH_PACKAGE", "PHG");
    sublist.set("EDGE_WEIGHT_DIM", "1");      // One weight per edge
    sublist.set("LB_APPROACH", "PARTITION");  // Build partition from scratch, in contrast to "REPARTITION"

    Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs = Teuchos::rcp(
        new Isorropia::Epetra::CostDescriber);

    costs->setGraphEdgeWeights(crs_ge_weights);

    Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner = Teuchos::rcp(
        new Isorropia::Epetra::Partitioner(initgraph_manip, costs, paramlist));

    Isorropia::Epetra::Redistributor rd(partitioner);
    Teuchos::RCP<Epetra_CrsGraph> bal_graph;

    //Use a try-catch block because Isorropia will throw an exception
    //if it encounters an error.

    try
    {
      bal_graph = rd.redistribute(*initgraph_manip);
    } catch (std::exception& exc)
    {
      std::cout << "Redistribute domain: Isorropia::Epetra::Redistributor threw "
          << "exception '" << exc.what() << std::endl;
      MPI_Finalize();
    }

    bal_graph->FillComplete();
    bal_graph->OptimizeStorage();

    return bal_graph;

}


/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> FSI::BlockMonolithic::SwitchDomains(
    Teuchos::RCP<Epetra_Map> rownodes, std::map<int, int>* nodeOwner,
    Teuchos::RCP<Epetra_CrsGraph> bal_graph, const Epetra_Comm& comm)
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
  comm.MaxAll(&nN,&maxNumOfNodes,1);

  int numMyFutureNodes = rownodes->NumMyElements();
  int myFutureNodeGIDs[maxNumOfNodes];    // contains GIDs of proc's nodes, later needed for map construction
  int* nodeIdsPtr = rownodes->MyGlobalElements();

  for (int k=0; k<numMyFutureNodes; ++k)  // Just insert the old node IDs. Maybe not every proc will change something here.
    myFutureNodeGIDs[k] = nodeIdsPtr[k];

  int countNodes[numproc];                // We use this array to determine to which proc another proc has to send its nodes.

  double nonmatch = 0;      // How many of my nodes don't want to be with me?
  double ratio;             // nonmatch is double because later this ratio is needed.

  double numInterfaceNodes=0;
  int globalInterfaceNodes = 0;

  int interfaceproc = 0;    // Does this proc own interface nodes?

  bool checkProcAgain;      // Always the current nodes of the proc are checked. If proc 0
                            // switches with proc 1, we have to check proc 0 again (such that
                            // the old nodes of proc 1 are checked). If we continue to proc 1,
                            // the old nodes of proc 0 would be checked again.

  int sendToProc = -1;
  int numInterfaceNodesCopy = 0;


  for (int proc=0; proc<numproc; ++proc){
      checkProcAgain = true;
      while (checkProcAgain){
        checkProcAgain = false;
        numInterfaceNodes = 0;
        nonmatch = 0;
        ratio = 0.0;
        for (int i=0; i<numproc; ++i){
          countNodes[i]=0;
        }

        int numNodes = numMyFutureNodes;
        int increase;
        int increaseCopy;
        int wishProc;         // the proc a node wants to belong to
        int wishProcCopy;
        int nodeIndex;

        comm.Broadcast(&numNodes,1,proc);
        for (int i=0; i<numNodes; ++i){
          if (myrank==proc){
            nodeIndex = myFutureNodeGIDs[i];
          }
          comm.Broadcast(&nodeIndex,1,proc);
          increase = -1;
          wishProc = -1;
          try {                             // Check if node is interface node. If yes:
            wishProc = nodeOwner->at(nodeIndex);     // Save procID to which node wants to belong.
            if (wishProc!=proc){
              increase = 1; }               // Node is NOT on the proc where it needs to be.
            else {
              increase = 0; }               // Node is on the proc where it belongs.
            }
          catch (std::exception& exc) {
          }
          for (int sender=0; sender<numproc; ++sender){
            increaseCopy = increase;
            wishProcCopy = wishProc;
            comm.Broadcast(&increaseCopy,1,sender);
            comm.Broadcast(&wishProcCopy,1,sender);
            if (increaseCopy != -1){
              if (myrank == proc){
                nonmatch += increaseCopy;
                ++countNodes[wishProcCopy]; // Save the info to which proc this nodes wants to belong.
                                            // If the node has to be sent away it will be sent to the proc with
                                            // the highest number of countNodes.
                ++numInterfaceNodes;
              }
            }
          }
        } //for (int i=0; i<numNodes; ++i){

        // Check if I have to give away my nodes and if yes, to which proc.

        if (myrank==proc && numInterfaceNodes!=0){
          ratio = nonmatch/numInterfaceNodes;
          int sendToProcCopy = sendToProc;
          if (ratio > 0.5) {  // I need to give my nodes to another processor.
            sendToProc = 0;
            for (int p=1; p<numproc; ++p){
              if (countNodes[p]>countNodes[sendToProc])
                sendToProc = p;
            }
            if (sendToProc == sendToProcCopy){                // Case: The proc wants to send back its nodes to
              if (numInterfaceNodes <= numInterfaceNodesCopy) // the proc from which it received them. Give them to the proc with
                sendToProc = -1;                              // the biggest share of the opposing patch. Explanation why we need this below.
            }
            else
              interfaceproc = 0;
          }
          else{
            sendToProc = -1;  // I can keep my nodes.
            interfaceproc = 1;
          }
        }
        else{
          sendToProc = -1;
        }

        /*
         * Explanation of the case: Proc wants to send back its nodes to the proc from which it received them
         *
         * Some problems are sensible to the selected edge weights in the graph which results in a bad partition
         * where maybe one large patch on the undistributed domain is opposed to two small patches which together
         * match the big domain. Without detecting this case this results in an infinite loop sending the nodes
         * from one proc to the other and back. Current implementation:
         * Check which proc has the larger share of the opposing patch. Still not a good distribution but the best
         * we can achieve with the given decomposition.
         */

        numInterfaceNodesCopy = numInterfaceNodes;            // Save this value in the copy-variable until next loop iteration.

//        if (comm_.MyPID()==proc)
//          std::cout<<"\nProc "<<proc<<" has to switch with "<<sendToProc<<", it is ratio = "<<ratio<<std::endl;

        comm.Broadcast(&sendToProc,1,proc);

        // Now do the actual switch between two procs, proc and sendToProc.

        int numNodesForwards;   // Number of nodes given from proc to sendToProc.
        int numNodesBackwards;  // Number of nodes given from sendToProc to proc.

        int nodeForwards;       // Nodes given from proc to sendToProc.
        int nodeBackwards;      // Nodes given from sendToProc to proc.

        if (sendToProc != -1 && sendToProc != proc){
          for (int i=0; i<maxNumOfNodes; ++i){
            if (myrank == proc){
              nodeForwards = myFutureNodeGIDs[i];
            }
            comm.Broadcast(&nodeForwards,1,proc);
            if (myrank == sendToProc){
              nodeBackwards = myFutureNodeGIDs[i];
              myFutureNodeGIDs[i] = nodeForwards;
            }
            comm.Broadcast(&nodeBackwards,1,sendToProc);
            if (myrank == proc)
              myFutureNodeGIDs[i] = nodeBackwards;
          }
          if (myrank == proc){
            numNodesForwards = numMyFutureNodes;
          }
          comm.Broadcast(&numNodesForwards,1,proc);
          if (myrank == sendToProc){
            numNodesBackwards = numMyFutureNodes;
            numMyFutureNodes = numNodesForwards;
            interfaceproc = 1;
          }
          comm.Broadcast(&numNodesBackwards,1,sendToProc);
          if (myrank == proc)
            numMyFutureNodes = numNodesBackwards;

          checkProcAgain = true;  // Procs switched nodes. Now let the same proc check its new nodes.
        }

        comm.Broadcast(&numInterfaceNodes,1,proc);
        globalInterfaceNodes += numInterfaceNodes;
      }
    } //for (int proc=0; proc<numproc; ++proc){

    // Which procs own interface nodes?
    for (int p=0; p<numproc; ++p){
      int interfaceproccopy = interfaceproc;
      comm.Broadcast(&interfaceproccopy,1,p);
      if (interfaceproccopy)
        interfaceprocs_.push_back(p);
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

    int switch_err = switched_bal_graph->Export(*bal_graph,exporter,Insert);

    if (switch_err)
      dserror("Switch of domains not successful, error code %d!",switch_err);

    return switched_bal_graph;
}

/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::InsertDeletedEdges(
    std::map<int, std::list<int> >* deletedEdges,
    Teuchos::RCP<Epetra_Map> switched_rownodes,
    Teuchos::RCP<Epetra_CrsGraph> switched_bal_graph)
{
  // Insert deleted edges
  std::map<int,std::list<int> >::iterator delEdgeIt;
  int row;
  std::list<int> col_list;
  std::list<int>::iterator col_list_it;
  for(delEdgeIt = deletedEdges->begin(); delEdgeIt != deletedEdges->end(); delEdgeIt++){
    row = delEdgeIt->first;
    if (switched_rownodes->LID(row) != -1){
      col_list = delEdgeIt->second;
      int col_length = col_list.size();
      int col_indices[col_length];
      int k=0;
      for (col_list_it = col_list.begin(); col_list_it != col_list.end(); col_list_it++){
        col_indices[k] = *col_list_it;
        ++k;
      }
      int success = switched_bal_graph->InsertGlobalIndices(row,col_length,col_indices);
      if (success != 0)
        dserror("\nInsert global indices failed, error code %d!", success);
    }
  }

  switched_bal_graph->FillComplete();
  switched_bal_graph->OptimizeStorage();
}

/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::FindNodeRelatedToDof(
    std::map<int, DRT::Node*>* nodes, int gdofid,
    Teuchos::RCP<DRT::Discretization> discretization, int* re)
{
  re[0] = -2;         // code: the node cannot be found on this proc
  bool breakout = false;
  std::vector<int> dofs;
  std::map<int, DRT::Node*>::iterator nodeiterator;

  for (nodeiterator = nodes->begin(); nodeiterator != nodes->end();
      nodeiterator++)
  {
    discretization->Dof((const DRT::Node *) nodeiterator->second, dofs);
    for (int i = 0; i < (int) dofs.size(); ++i)
    {
      if (dofs[i] == gdofid)
      {
        re[0] = nodeiterator->first;
        re[1] = nodeiterator->second->Owner();
        breakout = true;
        break;
      }
    }
    if (breakout == true)
      break;
    dofs.clear();
  }
}
