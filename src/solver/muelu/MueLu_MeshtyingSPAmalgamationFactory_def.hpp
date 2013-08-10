/*
 * MueLu_MeshtyingSPAmalgamationFactory_def.hpp
 *
 *  Created on: 18.03.2013
 *      Author: wiesner
 */

#ifndef MUELU_MESHTYINGSPAMALGAMATIONFACTORY_DEF_HPP_
#define MUELU_MESHTYINGSPAMALGAMATIONFACTORY_DEF_HPP_

#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q1_2013

#include "MueLu_MeshtyingSPAmalgamationFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>

#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_AmalgamationInfo.hpp>
#include <MueLu_Aggregates.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MeshtyingSPAmalgamationFactory()
{

}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<const ParameterList> MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A");
  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~MeshtyingSPAmalgamationFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
  /*currentLevel.DeclareInput("A", AFact_.get(), this);
  currentLevel.DeclareInput("Aggregates", aggregatesFact_.get(), this);
  currentLevel.DeclareInput("UnAmalgamationInfo", amalgFact_.get(), this);

  currentLevel.DeclareInput("SlaveDofMap", MueLu::NoFactory::get(),this);*/
  Input(currentLevel, "A"); // blocked A
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

  LocalOrdinal  fullblocksize = 1;   // block dim for fixed size blocks
  GlobalOrdinal offset = 0;          // global offset of dof gids
  LocalOrdinal blockid = -1;         // block id in strided map
  LocalOrdinal nStridedOffset = 0;   // DOF offset for strided block id "blockid" (default = 0)
  LocalOrdinal stridedblocksize = fullblocksize; // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)

  // 1) check for blocking/striding information
  if(A->IsView("stridedMaps") &&
     Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
    Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap()) == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
    fullblocksize = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize(); // TODO shorten code
    offset   = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getOffset();
    blockid  = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getStridedBlockId();
    if (blockid > -1) {
      std::vector<size_t> stridingInfo = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getStridingData();
      for(size_t j=0; j<Teuchos::as<size_t>(blockid); j++) {
        nStridedOffset += stridingInfo[j];
      }
      stridedblocksize = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
    } else {
      stridedblocksize = fullblocksize;
    }
    oldView = A->SwitchToView(oldView);
    GetOStream(Runtime1, 0) << "MeshtyingSPAmalgamationFactory::Build(): found fullblocksize=" << fullblocksize << " and stridedblocksize=" << stridedblocksize << " from strided maps. offset=" << offset << std::endl;
    /*std::cout << "fullblocksize: " << fullblocksize << std::endl;
      std::cout << "offset: " << offset << std::endl;
      std::cout << "blockid: " << blockid << std::endl;
      std::cout << "nStridedOffset: " << nStridedOffset << std::endl;
      std::cout << "stridedblocksize: " << stridedblocksize << std::endl;*/
  } else GetOStream(Warnings0, 0) << "MeshtyingSPAmalgamationFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;
  // TODO: maybe no striding information on coarser levels -> misuse nullspace vector?

  // 2) prepare maps for amalgamated graph of A and
  //    setup unamalgamation information

  RCP<std::vector<GlobalOrdinal> > gNodeIds; // contains global node ids on current proc
  gNodeIds = Teuchos::rcp(new std::vector<GlobalOrdinal>);
  gNodeIds->empty();

  // in nodegid2dofgids_ for each node on the current proc a vector of
  // the corresponding DOFs gids is stored.
  // The map contains all nodes the current proc has connections to (including
  // nodes that are stored on other procs when there are off-diagonal entries in A)
  nodegid2dofgids_ = Teuchos::rcp(new std::map<GlobalOrdinal,std::vector<GlobalOrdinal> >);

  // use row map (this is only working if A is zero matrix, otherwise we would have to use the column map)
  Teuchos::RCP<const Map> colMap = A->getDomainMap();
  GlobalOrdinal cnt_amalRows = 0; // counts number of nodes (rows in amalgamated matrix) on current proc
  LocalOrdinal nColEle = Teuchos::as<LocalOrdinal>(A->getDomainMap()->getNodeNumElements());
  for(LocalOrdinal i=0; i<nColEle;i++) {
    // get global DOF id
    GlobalOrdinal gDofId = colMap->getGlobalElement(i);

    // translate DOFGid to node id
#if defined( HAVE_Trilinos_Q3_2013)
    GlobalOrdinal gNodeId = AmalgamationFactory::DOFGid2NodeId(gDofId, fullblocksize, offset, 0 /*indexBase*/);
#elif defined(HAVE_Trilinos_Q2_2013)
    GlobalOrdinal gNodeId = AmalgamationFactory::DOFGid2NodeId(gDofId, A, fullblocksize, offset, 0 /*indexBase*/);
#else
    GlobalOrdinal gNodeId = AmalgamationFactory::DOFGid2NodeId(gDofId, A, fullblocksize, offset);
#endif

    // gblockid -> gDofId/lDofId
    if(nodegid2dofgids_->count(gNodeId) == 0) {

      // current column DOF gDofId belongs to a node that has not been added
      // to nodeid2dofgids_ yet. Do it now and add ALL DOFs of node gNodeId to
      // unamalgamation information.
      // Note: we use offset and fullblocksize, ie. information from strided maps indirectly
      std::vector<GlobalOrdinal> DOFs;

      DOFs.reserve(stridedblocksize);
      for(LocalOrdinal k=0; k<stridedblocksize; k++) {
        DOFs.push_back(offset + gNodeId*fullblocksize + nStridedOffset + k);
      }

      (*nodegid2dofgids_)[gNodeId] = DOFs;

      if(A->getRowMap()->isNodeGlobalElement(gDofId)) {
        gNodeIds->push_back(gNodeId);
        cnt_amalRows++; // new local block row in amalgamated matrix graph
      }
    }
  }

  // store (un)amalgamation information on current level
  RCP<AmalgamationInfo> amalgamationData = rcp(new AmalgamationInfo());
  amalgamationData->SetAmalgamationParams(nodegid2dofgids_);
  amalgamationData->SetNodeGIDVector(gNodeIds);
  amalgamationData->SetNumberOfNodes(cnt_amalRows);
  Set(currentLevel, "UnAmalgamationInfo", amalgamationData);

  //currentLevel.Set("Aggregates", aggregates, this);
}
} // namespace MueLu

#endif // Q1/2013
#endif // HAVE_MueLu



#endif /* MUELU_MESHTYINGSPAMALGAMATIONFACTORY_DEF_HPP_ */
