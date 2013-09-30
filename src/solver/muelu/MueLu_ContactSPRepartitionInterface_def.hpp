/*
 * MueLu_ContactSPRepartitionInterface_def.hpp
 *
 *  Created on: 11 Sep 2013
 *      Author: wiesner
 */

#ifndef MUELU_CONTACTSPREPARTITIONINTERFACE_DEF_HPP_
#define MUELU_CONTACTSPREPARTITIONINTERFACE_DEF_HPP_

#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q3_2013

#include "MueLu_ContactSPRepartitionInterface_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Utilities.hpp"


namespace MueLu {

 template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
 RCP<const ParameterList> ContactSPRepartitionInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("A",                    Teuchos::null, "Factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("AmalgamatedPartition", Teuchos::null, "(advanced) Factory generating the AmalgamatedPartition (e.g. an IsorropiaInterface)");
    //validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo",   Teuchos::null, "Generating factory of UnAmalgamationInfo");
    //validParamList->set< RCP<const FactoryBase> >("LagrNodeId2DispNodeId",Teuchos::null, "Generating factory for the mapping of Lagrange multiplier nodes to displacement DOF nodes. This has to be set by the user!");

    return validParamList;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactSPRepartitionInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "AmalgamatedPartition");
    //Input(currentLevel, "UnAmalgamationInfo");
    //Input(currentLevel, "LagrNodeId2DispNodeId");
    //currentLevel.DeclareInput("LagrNodeId2DispNodeId", MueLu::NoFactory::get() /*GetFactory("LagrNodeId2DispNodeId").get()*/,this);

  } //DeclareInput()

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ContactSPRepartitionInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &level) const {
    FactoryMonitor m(*this, "Build", level);
    level.print(GetOStream(Statistics0,0));
    // extract blocked operator A from current level
    RCP<Matrix> A = Get< RCP<Matrix> >     (level, "A");
    RCP<const Teuchos::Comm< int > > comm = A->getRowMap()->getComm();
    const int myRank = comm->getRank();

    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), false);
    ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);

    // fill decomposition vector
    for(LO i = 0; i<decomposition->getMap()->getNodeNumElements(); i++) {
      decompEntries[i] = myRank; // just a try -> put all to proc zero
    }

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    //decomposition->describe(*fos, Teuchos::VERB_EXTREME);

    Set(level, "Partition", decomposition);


#if 0 // not working
    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A);
    TEUCHOS_TEST_FOR_EXCEPTION(bA==Teuchos::null, Exceptions::BadCast, "MueLu::ContactSPRepartitionInterface::Build: input matrix A is not of type BlockedCrsMatrix! error.");

    // plausibility check
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows() != 2, Exceptions::RuntimeError, "MueLu::ContactSPRepartitionInterface::Build: number of block rows of A is not equal 2. error.");
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Cols() != 2, Exceptions::RuntimeError, "MueLu::ContactSPRepartitionInterface::Build: number of block cols of A is not equal 2. error.");

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));


    RCP<const Teuchos::Comm< int > > comm = A->getRowMap()->getComm();
    const int myRank = comm->getRank();

    //
    // ok, build
    //

    // pick out subblocks of block matrix
    RCP<CrsMatrix> A00 = bA->getMatrix(0,0);
    RCP<CrsMatrix> A01 = bA->getMatrix(0,1);
    RCP<CrsMatrix> A11 = bA->getMatrix(1,1);

    std::cout << "rowMap" << std::endl;
    A11->getRowMap()->describe(GetOStream(Statistics0,0),Teuchos::VERB_EXTREME);

    std::cout << "rangeMap" << std::endl;
    A11->getRangeMap()->describe(GetOStream(Statistics0,0),Teuchos::VERB_EXTREME);

    std::cout << "domainMap" << std::endl;
    A11->getDomainMap()->describe(GetOStream(Statistics0,0),Teuchos::VERB_EXTREME);

    std::cout << "colMap" << std::endl;
    A11->getColMap()->describe(GetOStream(Statistics0,0),Teuchos::VERB_EXTREME);

    // determine block information for displacement blocks
    // disp_offset usually is zero (default),
    // disp_blockdim is 2 or 3 (for 2d or 3d problems) on the finest level (# displacement dofs per node) and
    // disp_blockdim is 3 or 6 (for 2d or 3d problems) on coarser levels (# nullspace vectors)
    LocalOrdinal  disp_blockdim = 1;         // block dim for fixed size blocks
    GlobalOrdinal disp_offset   = 0;         // global offset of dof gids
    if(Teuchos::rcp_dynamic_cast<const StridedMap>(bA->getRangeMap(0)) != Teuchos::null) {
      RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(bA->getRangeMap(0));
      TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::ContactSPAggregationFactory::Build(): cast to strided row map failed.");
      disp_blockdim = strMap->getFixedBlockSize(); disp_offset   = strMap->getOffset();
      GetOStream(Statistics, 0) << "ContactSPRepartitionInterface::Build():" << " found disp_blockdim=" << disp_blockdim << " from strided maps. disp_offset=" << disp_offset << std::endl;
    }

    // determine block information for Lagrange multipliers
    // lagr_offset usually > zero (set by domainOffset for Ptent11Fact)
    // lagr_blockdim is disp_blockdim (for 2d or 3d problems) on the finest level (1 Lagrange multiplier per displacement dof) and
    // lagr_blockdim is 2 or 3 (for 2d or 3d problems) on coarser levels (same as on finest level, whereas there are 3 or 6 displacement dofs per node)
    LocalOrdinal  lagr_blockdim = disp_blockdim;         // block dim for fixed size blocks
    GlobalOrdinal lagr_offset   = 1000;         // global offset of dof gids (TODO fix this)
    if(Teuchos::rcp_dynamic_cast<const StridedMap>(bA->getRangeMap(1)) != Teuchos::null) {
      RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(bA->getRangeMap(1));
      TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::ContactSPAggregationFactory::Build(): cast to strided row map failed.");
      lagr_blockdim = strMap->getFixedBlockSize(); lagr_offset   = strMap->getOffset();
      GetOStream(Statistics, 0) << "ContactSPRepartitionInterface::Build():" << " found lagr_blockdim=" << lagr_blockdim << " from strided maps. lagr_offset=" << lagr_offset << std::endl;
    }

    // fetch map with slave Dofs from Level
    // slaveDofMap contains all global slave displacement DOF ids on the current level
    RCP<const Map> slaveDofMap = level.Get< RCP<const Map> >("SlaveDofMap",MueLu::NoFactory::get());

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

    for (size_t r = 0; r < slaveDofMap->getNodeNumElements(); r++) {  // todo improve me: loop only over node (skip all dofs in between..)
      // obtain global slave displacement DOF id
      GlobalOrdinal disp_grid = A01->getRowMap()->getGlobalElement(r);

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
        A01->getLocalRowView(r, lagr_indices, lagr_vals);

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

          TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<GlobalOrdinal>(local_lagrNodeId2dispNodeId.size())<lagr_nodeId-gMinLagrNodeId,Exceptions::BadCast,"MueLu::ContactSPRepartitionInterface::Build(): lagrNodeId2dispNodeId.size()<lagr_nodeId-gMinLagrNodeId. error.");
          if(lagrNodeId2dispNodeId[lagr_nodeId-gMinLagrNodeId] == -1)
            local_lagrNodeId2dispNodeId[lagr_nodeId-gMinLagrNodeId] = disp_nodeId;
        }

    }

    GlobalOrdinal lagrNodeId2dispNodeId_size = Teuchos::as<GlobalOrdinal>(local_lagrNodeId2dispNodeId.size());
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,lagrNodeId2dispNodeId_size,&local_lagrNodeId2dispNodeId[0],&lagrNodeId2dispNodeId[0]);

    for(int k=0; k<lagrNodeId2dispNodeId_size; k++) {
      std::cout << k << "->" << lagrNodeId2dispNodeId[k] << std::endl;
    }

    // now we have lagrNodeId2dispNodeId

    // fetch amalgamated repartitioning information from level
    // amalgamated partition is based on displacement node map
    RCP<Xpetra::Vector<GO, LO, GO, NO> > amalgPartition = Get< RCP<Xpetra::Vector<GO, LO, GO, NO> > >(level, "AmalgamatedPartition");
    ArrayRCP<GO> amalgPartitionData    = amalgPartition->getDataNonConst(0);
    RCP<const Map> dispNodeMap         = amalgPartition->getMap();

    TEUCHOS_TEST_FOR_EXCEPTION(dispNodeMap->getNodeNumElements()*disp_blockdim != A00->getRowMap()->getNodeNumElements(), Exceptions::RuntimeError, "MueLu::ContactSPAggregationFactory::Build(): Inconsistency between nodeMap and dofMap");

    // vector which stores final (unamalgamated) repartitioning
    // based on the Lagrange DOF map
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A11->getRowMap(), false);
    ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);

    // fill decomposition vector
    for(LO i = 0; i<decomposition->getMap()->getNodeNumElements(); i++) {
      GO gLagrDofId = decomposition->getMap()->getGlobalElement(i);

      GO gLagrNodeId = AmalgamationFactory::DOFGid2NodeId(gLagrDofId, lagr_blockdim, lagr_offset, 0/*indexBase*/);

      GO gDispNodeId = lagrNodeId2dispNodeId[gLagrNodeId];
      LO lDispNodeId = dispNodeMap->getLocalElement(gDispNodeId);

      std::cout << "gLagrDofId = " << gLagrDofId << " gLagrNodeId = " << gLagrNodeId << " gDispNodeId=" << gDispNodeId << " lDispNodeId=" << lDispNodeId << " amalgPartition.size()=" << amalgPartition->getMap()->getNodeNumElements() << std::endl;

      decompEntries[i] = myRank; // just a try -> put all to proc zero
      //if(lNodeId != -1)
      //  decompEntries[i] = amalgPartitionData[lNodeId];
      //for(size_t j=0; j<stridedblocksize/*DOFs.size()*/; j++) {
        // transform global DOF ids to local DOF ids using rowMap
        // note: The vector decomposition is based on rowMap
        //LO lDofId = rowMap->getLocalElement(DOFs[j]);     // -> i doubt that we need this!

        // put the same domain id to all DOFs of the same node
        //decompEntries[i*stridedblocksize + j] = amalgPartitionData[lNodeId];
        //decompEntries[lDofId] = amalgPartitionData[i];
      //}
    }


    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    decomposition->describe(*fos, Teuchos::VERB_EXTREME);

    Set(level, "Partition", decomposition);

#endif // end of "not working"

#if 0 // outdated code
    ArrayRCP<GO> amalgPartitionData = amalgPartition->getDataNonConst(0);

    RCP<const Map> rowMap        = A->getRowMap();
    RCP<const Map> nodeMap       = amalgPartition->getMap();

    // extract amalgamation information from matrix A
    LO blockdim = 1;                          // block dim for fixed size blocks
    GO indexBase = rowMap->getIndexBase();    // index base of maps
    GO offset    = 0;
    LO blockid          = -1;  // block id in strided map
    LO nStridedOffset   = 0;   // DOF offset for strided block id "blockid" (default = 0)
    LO stridedblocksize = blockdim; // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)

    // 1) check for blocking/striding information
    //    fill above variables
    if(A->IsView("stridedMaps") &&
       Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
      Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
      RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
      TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::RepartitionInterface::Build: cast to strided row map failed.");
      blockdim = strMap->getFixedBlockSize();
      offset   = strMap->getOffset();
      blockid  = strMap->getStridedBlockId();
      if (blockid > -1) {
        std::vector<size_t> stridingInfo = strMap->getStridingData();
        for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
          nStridedOffset += stridingInfo[j];
        stridedblocksize = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);

      } else {
        stridedblocksize = blockdim;
      }
      oldView = A->SwitchToView(oldView);
      GetOStream(Statistics0, -1) << "RepartitionInterface::Build():" << " found blockdim=" << blockdim << " from strided maps (blockid=" << blockid << ", strided block size=" << stridedblocksize << "). offset=" << offset << std::endl;
    } else GetOStream(Statistics0, -1) << "RepartitionInterface::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;

    // vector which stores final (unamalgamated) repartitioning
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(rowMap, false);
    ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);

    //std::cout << "nodeMap->getNumElements()=" << nodeMap->getNodeNumElements() << " stridedblocksize " << stridedblocksize << " rowMap->getNodeNumElements()=" << rowMap->getNodeNumElements() << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<int>(nodeMap->getNodeNumElements())*stridedblocksize != Teuchos::as<int>(rowMap->getNodeNumElements()), Exceptions::RuntimeError, "Inconsistency between nodeMap and dofMap");

    RCP<std::map<GO,std::vector<GO> > > nodegid2dofgids = amalgInfo->GetGlobalAmalgamationParams();

#if 0
    // fill decomposition vector
    for(LO i = 0; i<decomposition->getMap()->getNodeNumElements(); i++) {
      GO gDofId = decomposition->getMap()->getGlobalElement(i);

      GO gNodeId = AmalgamationFactory::DOFGid2NodeId(gDofId, blockdim, offset, indexBase);
      LO lNodeId = nodeMap->getLocalElement(gNodeId);

      std::cout << "gDofId = " << gDofId << " gNodeId = " << gNodeId << " lNodeId = " << lNodeId << " amalgPartition.size()=" << amalgPartition->getMap()->getNodeNumElements() << std::endl;

      if(lNodeId != -1)
        decompEntries[i] = amalgPartitionData[lNodeId];
      //for(size_t j=0; j<stridedblocksize/*DOFs.size()*/; j++) {
        // transform global DOF ids to local DOF ids using rowMap
        // note: The vector decomposition is based on rowMap
        //LO lDofId = rowMap->getLocalElement(DOFs[j]);     // -> i doubt that we need this!

        // put the same domain id to all DOFs of the same node
        //decompEntries[i*stridedblocksize + j] = amalgPartitionData[lNodeId];
        //decompEntries[lDofId] = amalgPartitionData[i];
      //}
    }
#else
    // fill vector with information about partitioning
    // TODO: we assume simple block maps here
    // TODO: adapt this to usage of nodegid2dofgids
    for(LO i = 0; i<nodeMap->getNodeNumElements(); i++) {
      // not fully sure about this. We're filling local ids in the decomposition vector with
      // the results stored in array. The decomposition vector is created using the rowMap of A

      // transform local node id to global node id.
      //GO gNodeId = nodeMap->getGlobalElement(i);

      // extract global DOF ids that belong to gNodeId
      /*std::vector<GlobalOrdinal> DOFs = (*nodegid2dofgids)[gNodeId];
      for(size_t j=0; j<stridedblocksize; j++) {
        decompEntries[i*stridedblocksize + j] = myRank;
      }*/
      for(size_t j=0; j<stridedblocksize/*DOFs.size()*/; j++) {
        // transform global DOF ids to local DOF ids using rowMap
        // note: The vector decomposition is based on rowMap
        //LO lDofId = rowMap->getLocalElement(DOFs[j]);     // -> i doubt that we need this!

        // put the same domain id to all DOFs of the same node
        decompEntries[i*stridedblocksize + j] = amalgPartitionData[i];
        //decompEntries[lDofId] = amalgPartitionData[i];
      }

    }
#endif

    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    decomposition->describe(*fos, Teuchos::VERB_EXTREME);

    Set(level, "Partition", decomposition);
#endif
  } //Build()



} //namespace MueLu

#endif
#endif

#endif /* MUELU_CONTACTSPREPARTITIONINTERFACE_DEF_HPP_ */
