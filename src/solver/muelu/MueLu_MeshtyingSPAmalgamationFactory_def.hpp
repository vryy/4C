/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu amalgamation factory for meshtying
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_MESHTYINGSPAMALGAMATIONFACTORY_DEF_HPP_
#define MUELU_MESHTYINGSPAMALGAMATIONFACTORY_DEF_HPP_

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

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal,
      Node>::MeshtyingSPAmalgamationFactory()
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const ParameterList>
  MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList(
      const ParameterList& paramList) const
  {
    Teuchos::RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set<Teuchos::RCP<const FactoryBase>>(
        "A", Teuchos::null, "Generating factory of the matrix A");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal,
      Node>::~MeshtyingSPAmalgamationFactory()
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(
      Level& currentLevel) const
  {
    /*currentLevel.DeclareInput("A", AFact_.get(), this);
    currentLevel.DeclareInput("Aggregates", aggregatesFact_.get(), this);
    currentLevel.DeclareInput("UnAmalgamationInfo", amalgFact_.get(), this);

    currentLevel.DeclareInput("SlaveDofMap", MueLu::NoFactory::get(),this);*/
    Input(currentLevel, "A");  // blocked A
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
      Level& currentLevel) const
  {
    FactoryMonitor m(*this, "Build", currentLevel);

    Teuchos::RCP<Matrix> A = Get<Teuchos::RCP<Matrix>>(currentLevel, "A");

    LocalOrdinal fullblocksize = 1;   // block dim for fixed size blocks
    GlobalOrdinal offset = 0;         // global offset of dof gids
    LocalOrdinal blockid = -1;        // block id in strided map
    LocalOrdinal nStridedOffset = 0;  // DOF offset for strided block id "blockid" (default = 0)
    LocalOrdinal stridedblocksize =
        fullblocksize;  // size of strided block id "blockid" (default = fullblocksize, only if
                        // blockid!=-1 stridedblocksize <= fullblocksize)

    // This is the cleaned code for Trilinos_Q1_2014

    // 1) check for blocking/striding information
    if (A->IsView("stridedMaps") &&
        Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null)
    {
      Xpetra::viewLabel_t oldView =
          A->SwitchToView("stridedMaps");  // note: "stridedMaps are always non-overlapping
                                           // (correspond to range and domain maps!)
      TEUCHOS_TEST_FOR_EXCEPTION(
          Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap()) == Teuchos::null,
          Exceptions::BadCast, "MueLu::CoalesceFactory::Build: cast to strided row map failed.");
      fullblocksize = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())
                          ->getFixedBlockSize();  // TODO shorten code
      offset = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getOffset();
      blockid = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getStridedBlockId();
      if (blockid > -1)
      {
        std::vector<size_t> stridingInfo =
            Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getStridingData();
        for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
        {
          nStridedOffset += stridingInfo[j];
        }
        stridedblocksize = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
      }
      else
      {
        stridedblocksize = fullblocksize;
      }
      oldView = A->SwitchToView(oldView);
      GetOStream(Runtime1, 0) << "MeshtyingSPAmalgamationFactory::Build(): found fullblocksize="
                              << fullblocksize << " and stridedblocksize=" << stridedblocksize
                              << " from strided maps. offset=" << offset << std::endl;
    }
    else
      GetOStream(Warnings0, 0) << "MeshtyingSPAmalgamationFactory::Build(): no striding "
                                  "information available. Use blockdim=1 with offset=0"
                               << std::endl;

    // 2) prepare maps for amalgamated graph of A and
    //    setup unamalgamation information

    Teuchos::RCP<std::vector<GlobalOrdinal>> gNodeIds;  // contains global node ids on current proc
    gNodeIds = Teuchos::rcp(new std::vector<GlobalOrdinal>);
    gNodeIds->empty();

    // use row map (this is only working if A is zero matrix, otherwise we would have to use the
    // column map)
    Teuchos::RCP<const Map> colMap = A->getDomainMap();
    LocalOrdinal nColEle = Teuchos::as<LocalOrdinal>(A->getDomainMap()->getNodeNumElements());
    for (LocalOrdinal i = 0; i < nColEle; i++)
    {
      // get global DOF id
      GlobalOrdinal gDofId = colMap->getGlobalElement(i);

      // translate DOFGid to node id
      GlobalOrdinal gNodeId =
          AmalgamationFactory::DOFGid2NodeId(gDofId, fullblocksize, offset, 0 /*indexBase*/);

      if (A->getRowMap()->isNodeGlobalElement(gDofId))
      {
        gNodeIds->push_back(gNodeId);
      }
    }

    // remove duplicates
    std::sort(gNodeIds->begin(), gNodeIds->end());
    gNodeIds->erase(std::unique(gNodeIds->begin(), gNodeIds->end()), gNodeIds->end());

    // store (un)amalgamation information on current level
    // TODO fix me
    Teuchos::RCP<AmalgamationInfo> amalgamationData =
        Teuchos::null;  // rcp(new AmalgamationInfo(gNodeIds,colMap,fullblocksize,
                        // offset,blockid,nStridedOffset,stridedblocksize));
    Set(currentLevel, "UnAmalgamationInfo", amalgamationData);
  }
}  // namespace MueLu

#endif /* MUELU_MESHTYINGSPAMALGAMATIONFACTORY_DEF_HPP_ */
