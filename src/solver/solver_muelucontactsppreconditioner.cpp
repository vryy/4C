/*
 * solver_muelucontactsppreconditioner.cpp
 *
 *  Created on: Sep 23, 2012
 *      Author: tobias
 */

#ifdef HAVE_MueLu

#include "../drt_lib/drt_dserror.H"

#include <MueLu_ConfigDefs.hpp>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
//#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_StridedEpetraMap.hpp>

// MueLu
#include <MueLu.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_IfpackSmoother.hpp>
#include <MueLu_SmootherPrototype_decl.hpp>

#include <MueLu_SubBlockAFactory.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
//#include <MueLu_UCAggregationFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_PgPFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_NullspaceFactory.hpp>
//#include <MueLu_SegregationAFilterFactory.hpp>
#include <MueLu_SegregationATransferFactory.hpp> // TODO remove me
#include <MueLu_Aggregates.hpp>
#include "MueLu_AggStatTransferFactory.hpp"
#include <MueLu_AggregationExportFactory.hpp>
#include <MueLu_BlockedPFactory.hpp>
#include <MueLu_BlockedGaussSeidelSmoother.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_SchurComplementFactory.hpp>
#include <MueLu_BraessSarazinSmoother.hpp>


#include <MueLu_MLParameterListInterpreter.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include <MueLu_UseShortNames.hpp>

#include <MueLu_EpetraOperator.hpp> // Aztec interface

#include "muelu/muelu_ContactAFilterFactory_decl.hpp"
#include "muelu/muelu_ContactTransferFactory_decl.hpp"
#include "muelu/muelu_ContactMapTransferFactory_decl.hpp"
#include "muelu/muelu_ContactASlaveDofFilterFactory_decl.hpp"
#include "muelu/muelu_ContactSPAggregationFactory_decl.hpp"
#include "muelu/MueLu_MyTrilinosSmoother_decl.hpp"

#include "solver_muelucontactsppreconditioner.H"

// BACI includes
#include "../linalg/linalg_blocksparsematrix.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuContactSpPreconditioner::MueLuContactSpPreconditioner( FILE * outfile, Teuchos::ParameterList & mllist )
  : PreconditionerType( outfile ),
    mllist_( mllist )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuContactSpPreconditioner::Setup( bool create,
                                              Epetra_Operator * matrix,
                                              Epetra_MultiVector * x,
                                              Epetra_MultiVector * b )
{
  std::cout << "call MueLuContactSpPreconditioner::Setup" << std::endl;

  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    // free old matrix first
    P_ = Teuchos::null;

    // adapt ML null space for contact/meshtying/constraint problems
    Teuchos::RCP<BlockSparseMatrixBase> A = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp( matrix, false ));
    if (A==Teuchos::null) dserror("matrix is not a BlockSparseMatrix");


    // prepare nullspace vector for MueLu (block A11 only)
    int numdf = mllist_.get<int>("PDE equations",-1);
    int dimns = mllist_.get<int>("null space: dimension",-1);
    if(dimns == -1 || numdf == -1) dserror("Error: PDE equations or null space dimension wrong.");

    // create a Teuchos::Comm from EpetraComm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::toXpetra(A->RangeMap(0).Comm());

    // get contact information
    Teuchos::RCP<Epetra_Map> epMasterDofMap = mllist_.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::MasterDofMap");
    Teuchos::RCP<Epetra_Map> epSlaveDofMap  = mllist_.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::SlaveDofMap");
    //Teuchos::RCP<Epetra_Map> epActiveDofMap = mllist_.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::ActiveDofMap");
    Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epSlaveDofMap  ));
    Teuchos::RCP<Xpetra::EpetraMap> xMasterDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epMasterDofMap  ));

    // create maps
    Teuchos::RCP<const Map> fullrangemap = Teuchos::rcp(new Xpetra::EpetraMap(Teuchos::rcpFromRef(A->FullRangeMap())));

    Teuchos::RCP<CrsMatrix> xA11 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0,0).EpetraMatrix()));
    Teuchos::RCP<CrsMatrix> xA12 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0,1).EpetraMatrix()));
    Teuchos::RCP<CrsMatrix> xA21 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1,0).EpetraMatrix()));
    Teuchos::RCP<CrsMatrix> xA22 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1,1).EpetraMatrix()));

    ///////////////////// DEBUG
    /*Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps> > t00 = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps>(xA11));
    Utils::Write("A00.dat", *t00);
    Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps> > t01 = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps>(xA12));
    Utils::Write("A01.dat", *t01);
    Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps> > t10 = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps>(xA21));
    Utils::Write("A10.dat", *t10);
    Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps> > t11 = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps>(xA22));
    Utils::Write("A11.dat", *t11);*/
    ///////////////////// DEBUG

    ///////////////////// EXPERIMENTAL
    std::vector<size_t> stridingInfo1;
    stridingInfo1.push_back(numdf);
    Teuchos::RCP<Xpetra::StridedEpetraMap> strMap1 = Teuchos::rcp(new Xpetra::StridedEpetraMap(Teuchos::rcpFromRef(A->Matrix(0,0).EpetraMatrix()->RowMap()), stridingInfo1, -1 /* stridedBlock */, 0 /*globalOffset*/));
    std::vector<size_t> stridingInfo2;
    stridingInfo2.push_back(numdf); // we have numdf Lagrange multipliers per node at the contact interface!
    Teuchos::RCP<Xpetra::StridedEpetraMap> strMap2 = Teuchos::rcp(new Xpetra::StridedEpetraMap(Teuchos::rcpFromRef(A->Matrix(1,1).EpetraMatrix()->RowMap()), stridingInfo2, -1 /* stridedBlock */, 0 /*globalOffset*/));
    /*std::cout << *strMap1 << std::endl;
    std::cout << *strMap2 << std::endl;*/
    ///////////////////// EXPERIMENTAL

    std::vector<Teuchos::RCP<const Map> > xmaps;
    xmaps.push_back(strMap1);
    xmaps.push_back(strMap2);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(fullrangemap,xmaps);

    // build blocked Xpetra operator
    Teuchos::RCP<BlockedCrsMatrix> bOp = Teuchos::rcp(new BlockedCrsMatrix(map_extractor,map_extractor,10));
    bOp->setMatrix(0,0,xA11);
    bOp->setMatrix(0,1,xA12);
    bOp->setMatrix(1,0,xA21);
    bOp->setMatrix(1,1,xA22);
    bOp->fillComplete();

    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap = xA11->getRowMap();

    Teuchos::RCP<MultiVector> nspVector11 = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
    Teuchos::RCP<std::vector<double> > nsdata = mllist_.get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

    for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
      Teuchos::ArrayRCP<Scalar> nspVector11i = nspVector11->getDataNonConst(i);
      const size_t myLength = nspVector11->getLocalLength();
      for(size_t j=0; j<myLength; j++) {
        nspVector11i[j] = (*nsdata)[i*myLength+j];
      }
    }

    ///////////////////////////////////////////////////////////////////////
    // special aggregation strategy
    //   - use 1pt aggregates for slave nodes
    ///////////////////////////////////////////////////////////////////////

    // number of node rows (only displacement dofs)
    const LocalOrdinal nDofRows = strMap1->getNodeNumElements();

    // prepare aggCoarseStat
    // TODO rebuild node-based map
    // still problematich for reparitioning
    Teuchos::ArrayRCP<unsigned int> aggStat;
    if(nDofRows > 0) aggStat = Teuchos::arcp<unsigned int>(nDofRows/numdf);
    for(LocalOrdinal i=0; i<nDofRows; ++i) {
      aggStat[i/numdf] = MueLu::NodeStats::READY;
      /*GlobalOrdinal grid = strMap1->getGlobalElement(i);
      if(xSlaveDofMap->isNodeGlobalElement(grid))
        aggStat[i/numdf] = MueLu::NodeStats::ONEPT;
      if(xMasterDofMap->isNodeGlobalElement(grid))
        aggStat[i/numdf] = MueLu::NodeStats::ONEPT;*/

    }

    // create Hierarchy
    Teuchos::RCP<Hierarchy> H = rcp(new Hierarchy());
    H->setDefaultVerbLevel(Teuchos::VERB_EXTREME);
    H->SetMaxCoarseSize(10); // TODO fix me
    H->GetLevel(0)->Set("A",Teuchos::rcp_dynamic_cast<Matrix>(bOp));
    H->GetLevel(0)->Set("Nullspace1",nspVector11);
    H->GetLevel(0)->Set("coarseAggStat",aggStat);
    H->GetLevel(0)->Set("MasterDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xMasterDofMap));  // set map with active dofs
    H->GetLevel(0)->Set("SlaveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xSlaveDofMap));  // set map with active dofs


    Teuchos::RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));
    Teuchos::RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 1));

    // set up block 11
    Teuchos::RCP<AmalgamationFactory> amalgFact11 = Teuchos::rcp(new AmalgamationFactory(A11Fact));
    amalgFact11->setDefaultVerbLevel(Teuchos::VERB_EXTREME);
    Teuchos::RCP<CoalesceDropFactory> dropFact11 = Teuchos::rcp(new CoalesceDropFactory(A11Fact,amalgFact11));
    dropFact11->setDefaultVerbLevel(Teuchos::VERB_EXTREME);
    Teuchos::RCP<UncoupledAggregationFactory> UCAggFact11 = Teuchos::rcp(new UncoupledAggregationFactory(dropFact11));
    UCAggFact11->SetMinNodesPerAggregate(9); // 9
    UCAggFact11->SetMaxNeighAlreadySelected(1);
    UCAggFact11->SetOrdering(MueLu::AggOptions::GRAPH);
    Teuchos::RCP<TentativePFactory> Ptent11Fact = Teuchos::rcp(new TentativePFactory(UCAggFact11,amalgFact11)); // check me
    Ptent11Fact->setStridingData(stridingInfo1);
    Teuchos::RCP<TentativePFactory> P11Fact = Ptent11Fact;
    //P11Fact->setStridedBlockId(0); // declare this P11Fact to be the transfer operator for the velocity dofs
    Teuchos::RCP<TransPFactory> R11Fact = Teuchos::rcp(new TransPFactory(P11Fact));
    Teuchos::RCP<NullspaceFactory> nspFact11 = Teuchos::rcp(new NullspaceFactory("Nullspace1",P11Fact));

    //////////////////////////////// define factory manager for (1,1) block
    Teuchos::RCP<FactoryManager> M11 = Teuchos::rcp(new FactoryManager());
    M11->SetFactory("A", A11Fact);
    M11->SetFactory("P", P11Fact);
    M11->SetFactory("Ptent", Ptent11Fact);
    M11->SetFactory("R", R11Fact);
    M11->SetFactory("Nullspace", nspFact11);
    M11->SetFactory("Ptent", P11Fact);
    M11->SetIgnoreUserData(true);               // always use data from factories defined in factory manager


    // create default nullspace
    int numPDEs = numdf;
    //GetOStream(MueLu::Runtime1, 0) << "Generating canonical nullspace: dimension = " << numPDEs << std::endl;
    Teuchos::RCP<MultiVector> nspVector22 = MultiVectorFactory::Build(xA22->getRowMap(), numPDEs);

    for (int i=0; i<numPDEs; ++i) {
      Teuchos::ArrayRCP<Scalar> nsValues22 = nspVector22->getDataNonConst(i);
      int numBlocks = nsValues22.size() / numPDEs;
      for (int j=0; j< numBlocks; ++j) {
        nsValues22[j*numPDEs + i] = 1.0;
      }
    }

    H->GetLevel(0)->Set("Nullspace2",nspVector22);

    // use TentativePFactory
    Teuchos::RCP<AmalgamationFactory> amalgFact22 = Teuchos::rcp(new AmalgamationFactory(A22Fact));
    Teuchos::RCP<MueLu::ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > UCAggFact22 = Teuchos::rcp(new MueLu::ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(UCAggFact11, amalgFact11));
    Teuchos::RCP<TentativePFactory> P22Fact = Teuchos::rcp(new TentativePFactory(UCAggFact22, amalgFact22));
    P22Fact->setStridingData(stridingInfo2);
    P22Fact->setDomainMapOffset(1000);
    //P22Fact->setStridedBlockId(0); // declare this P22Fact to be the transfer operator for the pressure dofs

    Teuchos::RCP<TransPFactory> R22Fact = Teuchos::rcp(new TransPFactory(P22Fact));

    Teuchos::RCP<NullspaceFactory> nspFact22 = Teuchos::rcp(new NullspaceFactory("Nullspace2",P22Fact));

    //////////////////////////////// define factory manager for (2,2) block
    Teuchos::RCP<FactoryManager> M22 = Teuchos::rcp(new FactoryManager());
    M22->SetFactory("A", A22Fact);
    M22->SetFactory("P", P22Fact);
    M22->SetFactory("R", R22Fact);
    M22->SetFactory("Aggregates", UCAggFact22);
    M22->SetFactory("Nullspace", nspFact22);
    M22->SetFactory("Ptent", P22Fact);
    M22->SetIgnoreUserData(true);               // always use data from factories defined in factory manager

    /////////////////////////////////////////// define blocked transfer ops
    Teuchos::RCP<BlockedPFactory> PFact = Teuchos::rcp(new BlockedPFactory(Teuchos::null)); // use row map index base from bOp
    PFact->AddFactoryManager(M11);
    PFact->AddFactoryManager(M22);

    Teuchos::RCP<GenericRFactory> RFact = Teuchos::rcp(new GenericRFactory(PFact));

    Teuchos::RCP<RAPFactory> AcFact = Teuchos::rcp(new RAPFactory(PFact, RFact));
    AcFact->SetRepairZeroDiagonal(true); // repair zero diagonal entries in Ac, that are resulting from Ptent with nullspacedim > ndofspernode

    // TODO add transfer factories!!!!
    Teuchos::RCP<MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cmTransFact3 = Teuchos::rcp(new MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", Ptent11Fact, MueLu::NoFactory::getRCP()));
    AcFact->AddTransferFactory(cmTransFact3);

    // create Braess-Sarazin smoother
    Scalar omega = 1.5;
    Teuchos::RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory(MueLu::NoFactory::getRCP(),omega));
    Teuchos::RCP<BraessSarazinSmoother> smootherPrototype = Teuchos::rcp(new BraessSarazinSmoother(1,omega)); // append SC smoother information

    Teuchos::RCP<SmootherFactory> smootherFact = Teuchos::rcp(new SmootherFactory(smootherPrototype));

    // SchurComplement smoother
    Teuchos::ParameterList SCList;
    /*SCList.set("relaxation: sweeps", (LO) 10);
    SCList.set("relaxation: damping factor", (SC) 0.6);
    SCList.set("relaxation: type", "Gauss-Seidel");
    Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::rcp(new TrilinosSmoother("RELAXATION",SCList,0,SFact));*/
    //Teuchos::RCP<SmootherPrototype> smoProtoSC = MueLu::GetIfpackSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("ILU", SCList,0,SFact);
    Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::rcp( new DirectSolver("Klu"/*"Umfpack"*/,Teuchos::ParameterList(),SFact) );
    Teuchos::RCP<SmootherFactory> SmooSCFact = Teuchos::rcp(new SmootherFactory(smoProtoSC));

    Teuchos::RCP<FactoryManager> MB = Teuchos::rcp(new FactoryManager());
    MB->SetFactory("A", SFact);
    MB->SetFactory("Smoother", SmooSCFact);
    MB->SetIgnoreUserData(true);
    smootherPrototype->SetFactoryManager(MB);  // add SC smoother information


    // main factory manager
    FactoryManager M;
    M.SetFactory("A", AcFact);
    M.SetFactory("P", PFact);
    M.SetFactory("R", RFact);
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("CoarseSolver", smootherFact);

    int maxLevels = 3;
    H->Setup(M,0,maxLevels);

    // set multigrid preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
  }

  //dserror("EXIT");

}



#endif // HAVE_MueLu

