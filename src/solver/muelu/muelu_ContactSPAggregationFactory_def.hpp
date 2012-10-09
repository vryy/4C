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

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ContactSPAggregationFactory(RCP<const FactoryBase> aggregatesFact, RCP<const FactoryBase> amalgFact)
: aggregatesFact_(aggregatesFact), amalgFact_(amalgFact), AFact_(MueLu::NoFactory::getRCP())
  {

  }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~ContactSPAggregationFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
  currentLevel.DeclareInput("A", AFact_.get(), this);
  currentLevel.DeclareInput("Aggregates", aggregatesFact_.get(), this);
  currentLevel.DeclareInput("UnAmalgamationInfo", amalgFact_.get(), this);
  currentLevel.DeclareInput("MasterDofMap", MueLu::NoFactory::get(),this); // TODO don't forget to transfer these maps!
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

  Monitor m(*this, "ContactSPAggregationFactory");


  RCP<Aggregates> aggs = currentLevel.Get<RCP<Aggregates> >("Aggregates", aggregatesFact_.get());
  RCP<LOVector> aggsvec = aggs->GetVertex2AggId();
  ArrayRCP< const LocalOrdinal > aggsdata = aggsvec->getData(0);
  std::cout << *aggsvec << std::endl;


  RCP<Matrix> Ain = currentLevel.Get< RCP<Matrix> >("A", AFact_.get());
  RCP<BlockedCrsMatrix> bOp = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ain);
  TEUCHOS_TEST_FOR_EXCEPTION(bOp==Teuchos::null, Exceptions::BadCast, "MueLu::ContactSPAggregationFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");

  const int myRank = bOp->getRangeMap()->getComm()->getRank();

  // pick out matrix block (0,1)
  RCP<CrsMatrix> A01 = bOp->getMatrix(0,1);

  // 1) check for blocking/striding information
  LocalOrdinal blockdim = 1;         // block dim for fixed size blocks
  GlobalOrdinal offset = 0;          // global offset of dof gids
  if(Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(0)) != Teuchos::null) {
    RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(bOp->getRangeMap(0));
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
    blockdim = strMap->getFixedBlockSize(); // TODO shorten code
    offset   = strMap->getOffset();
    GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build():" << " found blockdim=" << blockdim << " from strided maps. offset=" << offset << std::endl;
  } else GetOStream(Debug, 0) << "ContactSPAggregationFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;

  // loop over all rows
  // extract my values

  std::map<LocalOrdinal, std::vector<GlobalOrdinal> > laggId2gcids;
  std::map<LocalOrdinal, LocalOrdinal> dispAggId2lagAggId;
  std::map<LocalOrdinal, LocalOrdinal> nodeId2lagAggId;

  LocalOrdinal nLocalAggregates = 0;

  //fetch map with slave Dofs from Level
  RCP<const Map> slaveDofMap = currentLevel.Get< RCP<const Map> >("SlaveDofMap",MueLu::NoFactory::get());

  std::cout << "slaveDofMap " << std::endl;
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  slaveDofMap->describe(*fos,Teuchos::VERB_EXTREME);


  // handle slave dofs
  for (size_t r = 0; r < slaveDofMap->getNodeNumElements(); r++) {
    //std::cout << "r=" << r << std::endl;
    GlobalOrdinal grid = slaveDofMap->getGlobalElement(r);
    //std::cout << "grid("<<r<<")=" << grid << std::endl;

    if(A01->getRowMap()->isNodeGlobalElement(grid) ) {
      LocalOrdinal Alrid = A01->getRowMap()->getLocalElement(grid);

      // translate grid to nodeid
      GlobalOrdinal nodeId = AmalgamationFactory::DOFGid2NodeId(grid, Teuchos::null /* parameter not used */, blockdim, offset);
      LocalOrdinal lnodeId = aggsdata[Alrid/blockdim];

      Teuchos::ArrayView<const LocalOrdinal> indices;  // extract entries from local row
      Teuchos::ArrayView<const Scalar> vals;
      A01->getLocalRowView(Alrid, indices, vals);

      for(size_t i=0; i<(size_t)indices.size(); i++) {
        GlobalOrdinal gcid = A01->getColMap()->getGlobalElement(indices[i]);
        std::cout << "GRID: " << grid << " LID: " << Alrid << " node id: " << nodeId << " aggregate id: " << lnodeId << " GCID: " << gcid << " vals: " << vals[i] << std::endl;
        if(laggId2gcids.count(lnodeId) == 0) {
          std::vector<GlobalOrdinal> gcids;
          laggId2gcids[lnodeId] = gcids;
        }
        std::vector<GlobalOrdinal> & gcids = laggId2gcids[lnodeId];
        gcids.push_back(gcid);

        // displacementAggId2lagrangeMultAggId
        if(dispAggId2lagAggId.count(lnodeId) == 0) {
          dispAggId2lagAggId[lnodeId] = nLocalAggregates++;
        }
      }

      if(nodeId2lagAggId.count(nodeId) == 0) {
        nodeId2lagAggId[nodeId] = dispAggId2lagAggId[nodeId];
      }
    } //if A01->getRowMap()->isNodeGlobalElement()
  }

  typename std::map<LocalOrdinal, std::vector<GlobalOrdinal> >::iterator it;
  for(it = laggId2gcids.begin(); it != laggId2gcids.end(); it++) {
    std::cout << it->first << ": ";
    std::vector<GlobalOrdinal> & gcids = it->second;
    for(size_t t = 0; t<gcids.size(); t++) {
      std::cout << " " << gcids[t];
    }
    std::cout << std::endl;
  }

  // todo create nodemap for all slave nodes (usually non-overlapping) with
  // same distribution over processors than nodes

  // make sure that the unamalgamation information is correct (global node id -> global agg id)

  // generate "artificial nodes" for lagrange multipliers
  std::vector<GlobalOrdinal> lagMultNodes;
  for (size_t r = 0; r < bOp->getRangeMap(1)->getNodeNumElements(); r++) {
    //std::cout << "r=" << r << std::endl;
    GlobalOrdinal grid = bOp->getRangeMap(1)->getGlobalElement(r);
    //std::cout << "grid("<<r<<")=" << grid << std::endl;

    // translate grid to nodeid
    GlobalOrdinal nodeId = AmalgamationFactory::DOFGid2NodeId(grid, Teuchos::null /* parameter not used */, blockdim /*TODO */, offset);
    //std::cout << "nodeId: " << nodeId << std::endl;
    lagMultNodes.push_back(nodeId);
  }

  Teuchos::RCP<const Map > lagMultAggMap = MapFactory::Build(A01->getRowMap()->lib(),
                                                Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                lagMultNodes,
                                                /*indexBase*/ 0 /*TODO*/,
                                                A01->getRowMap()->getComm());

  std::cout << "lagMultAggMap " << std::endl;
  lagMultAggMap->describe(*fos,Teuchos::VERB_EXTREME);

  // build aggregates
  RCP<Aggregates> aggregates;

  // Build
  aggregates = rcp(new Aggregates(/*bOp->getRangeMap(1)*/lagMultAggMap));
  aggregates->setObjectLabel("UC (slave)");
  aggregates->SetNumAggregates(Teuchos::as<LocalOrdinal>(laggId2gcids.size())); // dont forget to set number of new aggregates

  // extract aggregate data structures to fill
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates->GetProcWinner()->getDataNonConst(0);

  // build new aggregates
  std::cout << "build " << laggId2gcids.size() << " more aggregates" << std::endl;

  for (size_t t = 0; t < slaveDofMap->getNodeNumElements(); t++) {
    GlobalOrdinal grid = slaveDofMap->getGlobalElement(t);
    //std::cout << "found grid " << grid << std::endl;
    if(A01->getRowMap()->isNodeGlobalElement(grid) ) {
      LocalOrdinal Alrid = A01->getRowMap()->getLocalElement(grid);

      // translate grid to nodeid
      GlobalOrdinal nodeId = AmalgamationFactory::DOFGid2NodeId(grid, Teuchos::null /* parameter not used */, blockdim, offset);
      LocalOrdinal laggId = aggsdata[Alrid/blockdim];

      std::cout << "found grid " << grid << " ~ Alrid: " << Alrid << " nodeId: " << nodeId << " laggId: " << laggId << std::endl;

      vertex2AggId[t] = dispAggId2lagAggId[laggId];
      //vertex2AggId[t] = nodeId2lagAggId[nodeId];
      procWinner[t] = myRank;
    }
  }


  // print results
  for(LocalOrdinal i = 0; i<vertex2AggId.size(); i++) {
    std::cout << "LID: " << i << " GID: " << aggregates->GetVertex2AggId()->getMap()->getGlobalElement(i) << " aggId: " << vertex2AggId[i] << " procWinner: " << procWinner[i] << std::endl;
  }

  currentLevel.Set("Aggregates", aggregates, this);
  aggregates->describe(GetOStream(Statistics0, 0), getVerbLevel());
}
} // namespace MueLu


#endif // HAVE_MueLu


#endif /* MUELU_CONTACTSPAGGREGATIONFACTORY_DEF_HPP_ */
