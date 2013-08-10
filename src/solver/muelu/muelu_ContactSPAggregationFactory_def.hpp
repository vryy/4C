/*
 * muelu_ContactSPAggregationFactory_def.hpp
 *
 *  Created on: Sep 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CONTACTSPAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_CONTACTSPAGGREGATIONFACTORY_DEF_HPP_

#ifdef HAVE_MueLu

#include "muelu_ContactSPAggregationFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>

#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_Aggregates.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ContactSPAggregationFactory(RCP<const FactoryBase> aggregatesFact, RCP<const FactoryBase> amalgFact)
: aggregatesFact_(aggregatesFact), amalgFact_(amalgFact), AFact_(MueLu::NoFactory::getRCP())
  {

  }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~ContactSPAggregationFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<const Teuchos::ParameterList> ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const Teuchos::ParameterList& paramList) const {
  Teuchos::RCP<Teuchos::ParameterList> validParamList = Teuchos::rcp(new Teuchos::ParameterList());

  // TODO remove aggregatesFact_, amalgFact_
  validParamList->set< Teuchos::RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set< Teuchos::RCP<const FactoryBase> >("Aggregates",     Teuchos::null, "Generating factory for aggregates");
  validParamList->set< Teuchos::RCP<const FactoryBase> >("UnAmalgamationInfo",Teuchos::null, "Generating factory for UnAmalgamationInfo.");
  validParamList->set< Teuchos::RCP<const FactoryBase> >("SlaveDofMap",MueLu::NoFactory::getRCP(), "Generating Factory for variable \"SlaveDofMap\"");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
  currentLevel.DeclareInput("A", AFact_.get(), this);
  currentLevel.DeclareInput("Aggregates", aggregatesFact_.get(), this);
  currentLevel.DeclareInput("UnAmalgamationInfo", amalgFact_.get(), this);
  //currentLevel.DeclareInput("MasterDofMap", MueLu::NoFactory::get(),this); // TODO don't forget to transfer these maps!
  currentLevel.DeclareInput("SlaveDofMap", MueLu::NoFactory::get(),this);

}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LOVector;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Matrix;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsMatrix;
  typedef Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> StridedMap;
  typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
  typedef MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Aggregates;
  //typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOOperator; //TODO
  //typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactoryClass;

  //Monitor m(*this, "ContactSPAggregationFactory");

  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  //slaveDofMap->describe(*fos,Teuchos::VERB_EXTREME);

#if 1
  // extract block matrix (must be a 2x2 block matrix)
  RCP<Matrix> Ain = currentLevel.Get< RCP<Matrix> >("A", AFact_.get());
  RCP<BlockedCrsMatrix> bOp = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ain);
  TEUCHOS_TEST_FOR_EXCEPTION(bOp==Teuchos::null, Exceptions::BadCast, "MueLu::ContactSPAggregationFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");

  // determine rank of current processor
  const Teuchos::RCP< const Teuchos::Comm< int > > comm = Ain->getRowMap()->getComm();
  const int myRank = comm->getRank();//bOp->getRangeMap()->getComm()->getRank();

  // pick out subblocks of block matrix
  RCP<CrsMatrix> A00 = bOp->getMatrix(0,0);
  RCP<CrsMatrix> A01 = bOp->getMatrix(0,1);

  // determine block information for displacement blocks
  // disp_offset usually is zero (default),
  // disp_blockdim is 2 or 3 (for 2d or 3d problems) on the finest level (# displacement dofs per node) and
  // disp_blockdim is 3 or 6 (for 2d or 3d problems) on coarser levels (# nullspace vectors)
  LocalOrdinal  disp_blockdim = 1;         // block dim for fixed size blocks
  GlobalOrdinal disp_offset   = 0;         // global offset of dof gids
  if(Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(0)) != Teuchos::null) {
    RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(0));
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::ContactSPAggregationFactory::Build(): cast to strided row map failed.");
    disp_blockdim = strMap->getFixedBlockSize(); disp_offset   = strMap->getOffset();
    //GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build():" << " found disp_blockdim=" << disp_blockdim << " from strided maps. disp_offset=" << disp_offset << std::endl;
  } //else GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build(): no striding information for displacement Dofs available. Use blockdim=1 with offset=0" << std::endl;

  // determine block information for Lagrange multipliers
  // lagr_offset usually > zero (set by domainOffset for Ptent11Fact)
  // lagr_blockdim is disp_blockdim (for 2d or 3d problems) on the finest level (1 Lagrange multiplier per displacement dof) and
  // lagr_blockdim is 2 or 3 (for 2d or 3d problems) on coarser levels (same as on finest level, whereas there are 3 or 6 displacement dofs per node)
  LocalOrdinal  lagr_blockdim = disp_blockdim;         // block dim for fixed size blocks
  GlobalOrdinal lagr_offset   = 1000;         // global offset of dof gids (TODO fix this)
  if(Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(1)) != Teuchos::null) {
    RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(1));
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::ContactSPAggregationFactory::Build(): cast to strided row map failed.");
    lagr_blockdim = strMap->getFixedBlockSize(); lagr_offset   = strMap->getOffset();
    //GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build():" << " found lagr_blockdim=" << lagr_blockdim << " from strided maps. lagr_offset=" << lagr_offset << std::endl;
  } else {
    //GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build(): no striding information for Lagrange multipliers available. Use lagr_blockdim=disp_blockdim=" << lagr_blockdim << " with lagr_offset=" << lagr_offset << std::endl;
  }

  // extract aggregates built using the displacement DOFs (from matrix block A00)
  RCP<Aggregates> dispAggs = currentLevel.Get<RCP<Aggregates> >("Aggregates", aggregatesFact_.get());
  RCP<LOVector> dispAggsVec = dispAggs->GetVertex2AggId();
  ArrayRCP< const LocalOrdinal > dispAggsData = dispAggsVec->getData(0);

  //dispAggsVec->describe(*fos,Teuchos::VERB_EXTREME);

  // fetch map with slave Dofs from Level
  // slaveDofMap contains all global slave displacement DOF ids on the current level
  RCP<const Map> slaveDofMap = currentLevel.Get< RCP<const Map> >("SlaveDofMap",MueLu::NoFactory::get());

  // generate global replicated mapping "lagrNodeId -> dispNodeId"
  RCP<const Map> lagrDofMap = A01->getDomainMap();
  GlobalOrdinal gMaxLagrNodeId = AmalgamationFactory::DOFGid2NodeId(
      lagrDofMap->getMaxAllGlobalIndex(),
#if defined(HAVE_Trilinos_Q3_2013)
#else
      Teuchos::null /* parameter not used */,
#endif
      lagr_blockdim,
      lagr_offset
#ifdef HAVE_Trilinos_Q2_2013
      , 0 /*indexBase*/
#endif
  );
  GlobalOrdinal gMinLagrNodeId = AmalgamationFactory::DOFGid2NodeId(
      lagrDofMap->getMinAllGlobalIndex(),
#if defined(HAVE_Trilinos_Q3_2013)
#else
      Teuchos::null /* parameter not used */,
#endif
      lagr_blockdim,
      lagr_offset
#ifdef HAVE_Trilinos_Q2_2013
      , 0 /*indexBase*/
#endif
  );

  // generate locally replicated vector for mapping Lagrange node ids to displacement node ids
  std::vector<GlobalOrdinal> lagrNodeId2dispNodeId(gMaxLagrNodeId-gMinLagrNodeId+1, -1);
  std::vector<GlobalOrdinal> local_lagrNodeId2dispNodeId(gMaxLagrNodeId-gMinLagrNodeId+1, -1);

  // generate locally replicated vector for mapping Lagrange node ids to displacement aggregate id.
  std::vector<GlobalOrdinal> lagrNodeId2dispAggId(gMaxLagrNodeId-gMinLagrNodeId+1, -1);
  std::vector<GlobalOrdinal> local_lagrNodeId2dispAggId(gMaxLagrNodeId-gMinLagrNodeId+1, -1);

  for (size_t r = 0; r < slaveDofMap->getNodeNumElements(); r++) {  // todo improve me: loop only over node (skip all dofs in between..)
    // obtain global slave displacement DOF id
    GlobalOrdinal disp_grid = slaveDofMap->getGlobalElement(r);

    if(A01->getRowMap()->isNodeGlobalElement(disp_grid) ) { // todo simplify me
      LocalOrdinal Alrid = A01->getRowMap()->getLocalElement(disp_grid);

      // get corresponding displacement aggregates.
      LocalOrdinal dispAggId = dispAggsData[Alrid/disp_blockdim];

      // translate displacement dof id to displacement node id
      GlobalOrdinal disp_nodeId = AmalgamationFactory::DOFGid2NodeId(
          disp_grid,
#if defined(HAVE_Trilinos_Q3_2013)
#else
          Teuchos::null /* parameter not used */,
#endif
          disp_blockdim,
          disp_offset
#ifdef HAVE_Trilinos_Q2_2013
          , 0 /*indexBase*/
#endif
          );

      Teuchos::ArrayView<const LocalOrdinal> lagr_indices;
      Teuchos::ArrayView<const Scalar> lagr_vals;
      A01->getLocalRowView(Alrid, lagr_indices, lagr_vals);

      for(size_t i = 0; i<Teuchos::as<size_t>(lagr_indices.size()); i++) {
        GlobalOrdinal lagr_gcid = A01->getColMap()->getGlobalElement(lagr_indices[i]);
        GlobalOrdinal lagr_nodeId = AmalgamationFactory::DOFGid2NodeId(
            lagr_gcid,
#if defined(HAVE_Trilinos_Q3_2013)
#else
            Teuchos::null /* parameter not used */,
#endif
            lagr_blockdim,
            lagr_offset
#ifdef HAVE_Trilinos_Q2_2013
            , 0 /*indexBase*/
#endif
            );

        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<GlobalOrdinal>(local_lagrNodeId2dispNodeId.size())<lagr_nodeId-gMinLagrNodeId,Exceptions::BadCast,"MueLu::ContactSPAggregationFactory::Build(): lagrNodeId2dispNodeId.size()<lagr_nodeId-gMinLagrNodeId. error.");
        if(lagrNodeId2dispNodeId[lagr_nodeId-gMinLagrNodeId] == -1)
          local_lagrNodeId2dispNodeId[lagr_nodeId-gMinLagrNodeId] = disp_nodeId;
        //else std::cout << "PROC: " << myRank << " lagr_nodeId " << lagr_nodeId << " is already connected to lagrange nodeId " <<  lagrNodeId2dispNodeId[lagr_nodeId-gMinLagrNodeId] << ". Ignore new dispNodeId: " << disp_nodeId << std::endl;

        local_lagrNodeId2dispAggId[lagr_nodeId-gMinLagrNodeId] = dispAggId; // todo: only if dispAggId belongs to current proc. these are local agg ids!
      }
    }
  }

  GlobalOrdinal lagrNodeId2dispNodeId_size = Teuchos::as<GlobalOrdinal>(local_lagrNodeId2dispNodeId.size());
  Teuchos::reduceAll(*comm /**A01->getColMap()->getComm()*/,Teuchos::REDUCE_MAX,lagrNodeId2dispNodeId_size,&local_lagrNodeId2dispNodeId[0],&lagrNodeId2dispNodeId[0]);
  Teuchos::reduceAll(*comm /**A01->getColMap()->getComm()*/,Teuchos::REDUCE_MAX,lagrNodeId2dispNodeId_size,&local_lagrNodeId2dispAggId[0],&lagrNodeId2dispAggId[0]);

  //Teuchos::reduceAll(*comm /**A01->getColMap()->getComm()*/,Teuchos::REDUCE_MAX,lagrNodeId2dispNodeId_size,&lagrNodeId2dispNodeId[0],&lagrNodeId2dispNodeId[0]);
  //Teuchos::reduceAll(*comm /**A01->getColMap()->getComm()*/,Teuchos::REDUCE_MAX,lagrNodeId2dispNodeId_size,&lagrNodeId2dispAggId[0],&lagrNodeId2dispAggId[0]);
  //Teuchos::reduceAll(*comm /**A01->getColMap()->getComm()*/,Teuchos::REDUCE_MAX,lagrNodeId2dispNodeId_size,&lagrNodeId2dispNodeId[0],Teuchos::outArg(&lagrNodeId2dispNodeId[0]));
  //Teuchos::reduceAll(*comm /**A01->getColMap()->getComm()*/,Teuchos::REDUCE_MAX,lagrNodeId2dispNodeId_size,&lagrNodeId2dispAggId[0],Teuchos::outArg(&lagrNodeId2dispAggId[0]));

  //Teuchos::outArg()
  /*for(size_t t = 0; t<lagrNodeId2dispNodeId.size(); t++) {
    std::cout << "PROC: " << myRank << " t=" << t << " lagrNodeId=" << t+gMinLagrNodeId << " corr dispNodeId=" << lagrNodeId2dispNodeId[t] << " dispAggId=" << lagrNodeId2dispAggId[t] << std::endl;
  }*/

  // build node map for Lagrange multipliers
  // generate "artificial nodes" for lagrange multipliers
  // the node map is also used for defining the Aggregates for the lagrange multipliers
  std::vector<GlobalOrdinal> lagr_Nodes;
  for (size_t r = 0; r < bOp->getRangeMap(1)->getNodeNumElements(); r++) {
    // determine global Lagrange multiplier row Dof
    // generate a node id using the grid, lagr_blockdim and lagr_offset // todo make sure, that nodeId is unique and does not interfer with the displacement nodes
    GlobalOrdinal lagr_grid = bOp->getRangeMap(1)->getGlobalElement(r);
    GlobalOrdinal lagr_nodeId =
        AmalgamationFactory::DOFGid2NodeId(
          lagr_grid,
#if defined(HAVE_Trilinos_Q3_2013)
#else
          Teuchos::null /* parameter not used */,
#endif
          lagr_blockdim, lagr_offset
#ifdef HAVE_Trilinos_Q2_2013
          , 0 /*indexBase*/
#endif
          );
    lagr_Nodes.push_back(lagr_nodeId);
  }

  // remove all duplicates
  lagr_Nodes.erase(std::unique(lagr_Nodes.begin(),lagr_Nodes.end()),lagr_Nodes.end());

  // define node map for Lagrange multipliers
  Teuchos::RCP<const Map > lagr_NodeMap = MapFactory::Build(A01->getRowMap()->lib(),
                                                Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), // TODO fix me
                                                lagr_Nodes,
                                                A01->getRowMap()->getIndexBase(),
                                                /*A01->getRowMap()->getComm()*/comm);

  //lagr_NodeMap->describe(*fos,Teuchos::VERB_EXTREME);

  // build processor local map: dispAggId2lagrAggId
  // generate new aggregegate ids if necessary (independent on each processor)

  // Build aggregates using the lagrange multiplier node map
  RCP<Aggregates> aggregates = Teuchos::rcp(new Aggregates(lagr_NodeMap));
  aggregates->setObjectLabel("UC (slave)");
  //aggregates->SetNumAggregates(Teuchos::as<LocalOrdinal>(dispAggId2lagAggId.size())); // dont forget to set number of new aggregates

  // extract aggregate data structures to fill
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates->GetProcWinner()->getDataNonConst(0);


  // loop over local lagrange multiplier node ids
  LocalOrdinal nLocalAggregates = 0;
  std::map<GlobalOrdinal,LocalOrdinal> dispAggId2localLagrAggId;
  for(size_t k = 0; k < lagr_NodeMap->getNodeNumElements(); k++) {
    GlobalOrdinal lagrNodeId = lagr_NodeMap->getGlobalElement(k);
    GlobalOrdinal dispAggId  = lagrNodeId2dispAggId[lagrNodeId-gMinLagrNodeId];
    if(dispAggId2localLagrAggId.count(dispAggId) == 0)
      dispAggId2localLagrAggId[dispAggId] = nLocalAggregates++;
    vertex2AggId[k] = dispAggId2localLagrAggId[dispAggId];
    procWinner[k] = myRank;
  }

  aggregates->SetNumAggregates(nLocalAggregates);

  // TODO think about this. probably not necessary...
  //aggregates->AggregatesCrossProcessors(true);  // force TentativePFactory to build a column map

  //aggregates->GetVertex2AggId()->describe(*fos,Teuchos::VERB_EXTREME);
  //aggregates->GetProcWinner()->describe(*fos,Teuchos::VERB_EXTREME);

  currentLevel.Set("Aggregates", aggregates, this);

#else // outdated
  // extract block matrix (must be a 2x2 block matrix)
  RCP<Matrix> Ain = currentLevel.Get< RCP<Matrix> >("A", AFact_.get());
  RCP<BlockedCrsMatrix> bOp = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ain);
  TEUCHOS_TEST_FOR_EXCEPTION(bOp==Teuchos::null, Exceptions::BadCast, "MueLu::ContactSPAggregationFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");

  // determine rank of current processor
  const int myRank = bOp->getRangeMap()->getComm()->getRank();

  std::cout << "PROC " << myRank << " ContactSPAggregationFactory::Build 1" << std::endl;

  // pick out subblocks of block matrix
  RCP<CrsMatrix> A00 = bOp->getMatrix(0,0);
  RCP<CrsMatrix> A01 = bOp->getMatrix(0,1);

  // determine block information for displacement blocks
  // disp_offset usually is zero (default),
  // disp_blockdim is 2 or 3 (for 2d or 3d problems) on the finest level (# displacement dofs per node) and
  // disp_blockdim is 3 or 6 (for 2d or 3d problems) on coarser levels (# nullspace vectors)
  LocalOrdinal  disp_blockdim = 1;         // block dim for fixed size blocks
  GlobalOrdinal disp_offset   = 0;         // global offset of dof gids
  if(Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(0)) != Teuchos::null) {
    RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(0));
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::ContactSPAggregationFactory::Build(): cast to strided row map failed.");
    disp_blockdim = strMap->getFixedBlockSize(); disp_offset   = strMap->getOffset();
    GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build():" << " found disp_blockdim=" << disp_blockdim << " from strided maps. disp_offset=" << disp_offset << std::endl;
  } else GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build(): no striding information for displacement Dofs available. Use blockdim=1 with offset=0" << std::endl;

  // determine block information for Lagrange multipliers
  // lagr_offset usually > zero (set by domainOffset for Ptent11Fact)
  // lagr_blockdim is disp_blockdim (for 2d or 3d problems) on the finest level (1 Lagrange multiplier per displacement dof) and
  // lagr_blockdim is 2 or 3 (for 2d or 3d problems) on coarser levels (same as on finest level, whereas there are 3 or 6 displacement dofs per node)
  LocalOrdinal  lagr_blockdim = disp_blockdim;         // block dim for fixed size blocks
  GlobalOrdinal lagr_offset   = 1000;         // global offset of dof gids (TODO fix this)
  if(Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(1)) != Teuchos::null) {
    RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(1));
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::ContactSPAggregationFactory::Build(): cast to strided row map failed.");
    lagr_blockdim = strMap->getFixedBlockSize(); lagr_offset   = strMap->getOffset();
    GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build():" << " found lagr_blockdim=" << lagr_blockdim << " from strided maps. lagr_offset=" << lagr_offset << std::endl;
  } else {
    GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build(): no striding information for Lagrange multipliers available. Use lagr_blockdim=disp_blockdim=" << lagr_blockdim << " with lagr_offset=" << lagr_offset << std::endl;
  }

  std::cout << "PROC " << myRank << " ContactSPAggregationFactory::Build 2" << std::endl;

  // extract aggregates built using the displacement DOFs (from matrix block A00)
  RCP<Aggregates> dispAggs = currentLevel.Get<RCP<Aggregates> >("Aggregates", aggregatesFact_.get());
  RCP<LOVector> dispAggsVec = dispAggs->GetVertex2AggId();
  ArrayRCP< const LocalOrdinal > dispAggsData = dispAggsVec->getData(0);

  dispAggsVec->describe(*fos,Teuchos::VERB_EXTREME);

  // fetch map with slave Dofs from Level
  // slaveDofMap contains all global slave displacement DOF ids on the current level
  // we extract the corresponding aggregate ids of the slave nodes
  RCP<const Map> slaveDofMap = currentLevel.Get< RCP<const Map> >("SlaveDofMap",MueLu::NoFactory::get());

  // generate map displacement aggregate id -> new Lagrange multiplier aggregate id
  // count number of new aggregates for Lagrange multipliers that are to be built on the current proc
  LocalOrdinal nLocalAggregates = 0;  // number of new local aggregates for Lagrange multipliers
  std::map<LocalOrdinal, LocalOrdinal> dispAggId2lagAggId;
  for (size_t r = 0; r < slaveDofMap->getNodeNumElements(); r++) {
    // obtain global slave displacement DOF id
    GlobalOrdinal grid = slaveDofMap->getGlobalElement(r);

    if(A01->getRowMap()->isNodeGlobalElement(grid) ) { // todo simplify me
      LocalOrdinal Alrid = A01->getRowMap()->getLocalElement(grid);

      // translate grid to nodeid
      //GlobalOrdinal nodeId = AmalgamationFactory::DOFGid2NodeId(grid, Teuchos::null /* parameter not used */, blockdim, offset);
      LocalOrdinal dispAggId = dispAggsData[Alrid/disp_blockdim];

      // displacementAggId2lagrangeMultAggId
      if(dispAggId2lagAggId.count(dispAggId) == 0) {
        dispAggId2lagAggId[dispAggId] = nLocalAggregates++;  // a new aggregate has to be built for this Lagrange multiplier
        std::cout << "dispAggId: " << dispAggId << " dispAggId2lagAggId[dispAggId]: " << dispAggId2lagAggId[dispAggId] << std::endl;
      }
    } //if A01->getRowMap()->isNodeGlobalElement()
  }

  std::cout << "PROC " << myRank << " ContactSPAggregationFactory::Build 3" << std::endl;

  // loop over all local slave dofs
  // disp_node -> lagr_dofs
  std::map<LocalOrdinal, std::vector<GlobalOrdinal> > dispNodeId2lagNodeIds;
  for (size_t r = 0; r < slaveDofMap->getNodeNumElements(); r++) {
    // obtain global slave displacement DOF id
    GlobalOrdinal disp_grid = slaveDofMap->getGlobalElement(r);

    if(A01->getRowMap()->isNodeGlobalElement(disp_grid) ) { // todo simplify me
      LocalOrdinal Alrid = A01->getRowMap()->getLocalElement(disp_grid);

      // translate displacement dof id to displacement node id
      GlobalOrdinal disp_nodeId = AmalgamationFactory::DOFGid2NodeId(
          disp_grid, Teuchos::null /* parameter not used */,
          disp_blockdim,
          disp_offset);

      Teuchos::ArrayView<const LocalOrdinal> lagr_indices;
      Teuchos::ArrayView<const Scalar> lagr_vals;
      A01->getLocalRowView(Alrid, lagr_indices, lagr_vals);

      std::cout << "disp_nodeId=" << disp_nodeId << " lagNodeIds:";

      std::vector<GlobalOrdinal> lagr_nodeIds;
      lagr_nodeIds.empty();
      for(size_t i = 0; i<lagr_indices.size(); i++) {
        GlobalOrdinal lagr_gcid = A01->getColMap()->getGlobalElement(lagr_indices[i]);
        GlobalOrdinal lagr_nodeId = AmalgamationFactory::DOFGid2NodeId(
            lagr_gcid, Teuchos::null /* parameter not used */,
            lagr_blockdim,
            lagr_offset);
        lagr_nodeIds.push_back(lagr_nodeId);
        std::cout << lagr_nodeId << "(" << lagr_gcid << ")" << " ";
      }
      dispNodeId2lagNodeIds[disp_nodeId] = lagr_nodeIds;
      std::cout << std::endl;

    }
  }

  std::cout << "PROC " << myRank << " ContactSPAggregationFactory::Build 4" << std::endl;

  bOp->getRangeMap(1)->describe(*fos,Teuchos::VERB_EXTREME);

  // generate "artificial nodes" for lagrange multipliers
  std::vector<GlobalOrdinal> lagr_Nodes;
  for (size_t r = 0; r < bOp->getRangeMap(1)->getNodeNumElements(); r++) {
    // determine global Lagrange multiplier row Dof
    // generate a node id using the grid, lagr_blockdim and lagr_offset // todo make sure, that nodeId is unique and does not interfer with the displacement nodes
    GlobalOrdinal lagr_grid = bOp->getRangeMap(1)->getGlobalElement(r);
    GlobalOrdinal lagr_nodeId =
        AmalgamationFactory::DOFGid2NodeId(
          lagr_grid, Teuchos::null /* parameter not used */,
          lagr_blockdim, lagr_offset);
    lagr_Nodes.push_back(lagr_nodeId);
  }

  // remove all duplicates
  lagr_Nodes.erase(std::unique(lagr_Nodes.begin(),lagr_Nodes.end()),lagr_Nodes.end());

  std::cout << "PROC " << myRank << " ContactSPAggregationFactory::Build 5" << std::endl;

  // define node map for Lagrange multipliers
  size_t numGlobalElements = 0;
  sumAll(A01->getRowMap()->getComm(), lagr_Nodes.size(), numGlobalElements);
  Teuchos::RCP<const Map > lagr_NodeMap = MapFactory::Build(A01->getRowMap()->lib(),
                                                //Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), // TODO fix me
                                                numGlobalElements,
                                                lagr_Nodes,
                                                A01->getRowMap()->getIndexBase(),
                                                A01->getRowMap()->getComm());

  lagr_NodeMap->describe(*fos,Teuchos::VERB_EXTREME);

  std::cout << "PROC " << myRank << " ContactSPAggregationFactory::Build 6" << std::endl;

  // Build aggregates using the lagrange multiplier node map
  RCP<Aggregates> aggregates = Teuchos::rcp(new Aggregates(lagr_NodeMap));
  aggregates->setObjectLabel("UC (slave)");
  aggregates->SetNumAggregates(Teuchos::as<LocalOrdinal>(dispAggId2lagAggId.size())); // dont forget to set number of new aggregates

  // extract aggregate data structures to fill
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates->GetProcWinner()->getDataNonConst(0);

  std::cout << "PROC " << myRank << " ContactSPAggregationFactory::Build 7" << std::endl;

  slaveDofMap->describe(*fos,Teuchos::VERB_EXTREME);

  for (size_t r = 0; r < slaveDofMap->getNodeNumElements(); r++) {
    // obtain global slave displacement DOF id
    GlobalOrdinal disp_grid = slaveDofMap->getGlobalElement(r);

    if(A01->getRowMap()->isNodeGlobalElement(disp_grid) ) { // todo simplify me
      LocalOrdinal Alrid = A01->getRowMap()->getLocalElement(disp_grid);

      LocalOrdinal dispAggId = dispAggsData[Alrid/disp_blockdim];

      std::cout << "myRank: " << myRank << " r=" << r << " disp_grid=" << disp_grid << " Alrid=" << Alrid << " dispAggId=" << dispAggId << std::endl;

      GlobalOrdinal disp_NodeId = AmalgamationFactory::DOFGid2NodeId(
          disp_grid, Teuchos::null /* parameter not used */,
          disp_blockdim,
          disp_offset);

      std::cout << "myRank: " << myRank << " r=" << r << " disp_grid=" << disp_grid << " Alrid=" << Alrid << " dispAggId=" << dispAggId << " disp_NodeId=" << disp_NodeId << std::endl;

      std::vector<GlobalOrdinal> lagr_nodeIds = dispNodeId2lagNodeIds[disp_NodeId];

      std::cout << "myRank: " << myRank << " r=" << r << " disp_grid=" << disp_grid << " Alrid=" << Alrid << " dispAggId=" << dispAggId << " disp_NodeId=" << disp_NodeId << " lagr_nodeIds.size()=" << lagr_nodeIds.size() << std::endl;

      for(size_t i=0; i< lagr_nodeIds.size(); i++) {
        GlobalOrdinal lagr_NodeId = lagr_nodeIds[i];
        std::cout << "myRank: " << myRank << " r=" << r << " disp_grid=" << disp_grid << " Alrid=" << Alrid << " dispAggId=" << dispAggId << " disp_NodeId=" << disp_NodeId << " lagr_nodeId=" << lagr_NodeId << std::endl;
        LocalOrdinal lagr_NodeLID = lagr_NodeMap->getLocalElement(lagr_NodeId); // here it fails
        std::cout << "myRank: " << myRank << " r=" << r << " disp_grid=" << disp_grid << " Alrid=" << Alrid << " dispAggId=" << dispAggId << " disp_NodeId=" << disp_NodeId << " lagr_nodeLID=" << lagr_NodeLID << std::endl;
        vertex2AggId[lagr_NodeLID] = dispAggId2lagAggId[dispAggId]; // aggregate
        std::cout << "myRank: " << myRank << "lagr_nodeLID: " << lagr_NodeLID << " lagr_NodeId: " << lagr_NodeId << " dispAggId: " << dispAggId << " lagrAggId: " <<  dispAggId2lagAggId[dispAggId] << std::endl;
        procWinner[lagr_NodeLID] = myRank;
      }
    }
  }

  std::cout << "PROC " << myRank << " ContactSPAggregationFactory::Build 8" << std::endl;

  std::cout << "lagAggregates" << std::endl;
  aggregates->GetVertex2AggId()->describe(*fos,Teuchos::VERB_EXTREME);

  currentLevel.Set("Aggregates", aggregates, this);
  aggregates->describe(GetOStream(Statistics0, 0), getVerbLevel());

  std::cout << "PROC " << myRank << " ContactSPAggregationFactory::Build 9" << std::endl;


#endif

}
} // namespace MueLu


#endif // HAVE_MueLu


#endif /* MUELU_CONTACTSPAGGREGATIONFACTORY_DEF_HPP_ */
