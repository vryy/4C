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

#include <MueLu_SegregationATransferFactory.hpp> // TODO remove me
#include <MueLu_Aggregates.hpp>
#include "MueLu_AggStatTransferFactory.hpp"
#include <MueLu_AggregationExportFactory.hpp>
#include <MueLu_BlockedPFactory.hpp>
//#include <MueLu_BlockedGaussSeidelSmoother.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_SchurComplementFactory.hpp>
#include <MueLu_BraessSarazinSmoother.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_BlockedCoarseMapFactory.hpp>

#include <MueLu_MLParameterListInterpreter.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include <MueLu_UseShortNames.hpp>

#include <MueLu_EpetraOperator.hpp> // Aztec interface

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

    ///////////////////////////////////////////////////////
    // interpret ML parameters
    int maxLevels = 3;
    int verbosityLevel = 10;
    int maxCoarseSize = 100;
    if(mllist_.isParameter("max levels")) maxLevels = mllist_.get<int>("max levels");
    if(mllist_.isParameter("ML output"))  verbosityLevel = mllist_.get<int>("ML output");
    if(mllist_.isParameter("coarse: max size")) maxCoarseSize = mllist_.get<int>("coarse: max size");

    // translate verbosity parameter
    Teuchos::EVerbosityLevel eVerbLevel = Teuchos::VERB_NONE;
    if(verbosityLevel == 0)  eVerbLevel = Teuchos::VERB_NONE;
    if(verbosityLevel > 0 )  eVerbLevel = Teuchos::VERB_LOW;
    if(verbosityLevel > 4 )  eVerbLevel = Teuchos::VERB_MEDIUM;
    if(verbosityLevel > 7 )  eVerbLevel = Teuchos::VERB_HIGH;
    if(verbosityLevel > 9 )  eVerbLevel = Teuchos::VERB_EXTREME;

    ////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    // prepare nullspace vector for MueLu (block A11 only)
    ///////////////////////////////////////////////////////////////////////
    int numdf = mllist_.get<int>("PDE equations",-1);
    int dimns = mllist_.get<int>("null space: dimension",-1);
    if(dimns == -1 || numdf == -1) dserror("Error: PDE equations or null space dimension wrong.");

    ///////////////////////////////////////////////////////////////////////
    // create a Teuchos::Comm from EpetraComm
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::toXpetra(A->RangeMap(0).Comm());

    ///////////////////////////////////////////////////////////////////////
    // get contact information
    ///////////////////////////////////////////////////////////////////////
    //Teuchos::RCP<Epetra_Map> epMasterDofMap = mllist_.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::MasterDofMap");
    Teuchos::RCP<Epetra_Map> epSlaveDofMap  = mllist_.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::SlaveDofMap");
    //Teuchos::RCP<Epetra_Map> epActiveDofMap = mllist_.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::ActiveDofMap");
    Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epSlaveDofMap  ));
    //Teuchos::RCP<Xpetra::EpetraMap> xMasterDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epMasterDofMap  ));

    ///////////////////////////////////////////////////////////////////////
    // prepare maps for blocked operator
    ///////////////////////////////////////////////////////////////////////

    // create maps
    Teuchos::RCP<const Map> fullrangemap = Teuchos::rcp(new Xpetra::EpetraMap(Teuchos::rcpFromRef(A->FullRangeMap())));

    Teuchos::RCP<CrsMatrix> xA11 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0,0).EpetraMatrix()));
    Teuchos::RCP<CrsMatrix> xA12 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0,1).EpetraMatrix()));
    Teuchos::RCP<CrsMatrix> xA21 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1,0).EpetraMatrix()));
    Teuchos::RCP<CrsMatrix> xA22 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1,1).EpetraMatrix()));

    // define strided maps
    std::vector<size_t> stridingInfo1;
    stridingInfo1.push_back(numdf);
    Teuchos::RCP<Xpetra::StridedEpetraMap> strMap1 = Teuchos::rcp(new Xpetra::StridedEpetraMap(Teuchos::rcpFromRef(A->Matrix(0,0).EpetraMatrix()->RowMap()), stridingInfo1, -1 /* stridedBlock */, 0 /*globalOffset*/));
    std::vector<size_t> stridingInfo2;
    stridingInfo2.push_back(numdf); // we have numdf Lagrange multipliers per node at the contact interface!
    Teuchos::RCP<Xpetra::StridedEpetraMap> strMap2 = Teuchos::rcp(new Xpetra::StridedEpetraMap(Teuchos::rcpFromRef(A->Matrix(1,1).EpetraMatrix()->RowMap()), stridingInfo2, -1 /* stridedBlock */, 0 /*globalOffset*/));

    // build map extractor
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

    ///////////////////////////////////////////////////////////////////////
    // prepare nullspace for first block
    ///////////////////////////////////////////////////////////////////////

    // extract nullspace information from ML list
    Teuchos::RCP<MultiVector> nspVector11 = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xA11->getRowMap(),dimns,true);
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

    ///////////////////////////////////////////////////////////////////////
    // create Hierarchy
    ///////////////////////////////////////////////////////////////////////

    Teuchos::RCP<Hierarchy> H = Teuchos::rcp(new Hierarchy());
    H->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(eVerbLevel));
    H->SetMaxCoarseSize(Teuchos::as<Xpetra::global_size_t>(maxCoarseSize));
    H->GetLevel(0)->Set("A",Teuchos::rcp_dynamic_cast<Matrix>(bOp));
    H->GetLevel(0)->Set("Nullspace1",nspVector11);
    H->GetLevel(0)->Set("coarseAggStat",aggStat);
    //H->GetLevel(0)->Set("MasterDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xMasterDofMap));  // set map with active dofs
    H->GetLevel(0)->Set("SlaveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xSlaveDofMap));  // set map with active dofs


    Teuchos::RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));
    Teuchos::RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 1));

    ///////////////////////////////////////////////////////////////////////
    // set up block 11
    ///////////////////////////////////////////////////////////////////////

    Teuchos::RCP<AmalgamationFactory> amalgFact11 = Teuchos::rcp(new AmalgamationFactory(/*A11Fact*/));
    //amalgFact11->setDefaultVerbLevel(Teuchos::VERB_EXTREME);
    Teuchos::RCP<CoalesceDropFactory> dropFact11 = Teuchos::rcp(new CoalesceDropFactory(/*A11Fact,amalgFact11*/));
    //dropFact11->setDefaultVerbLevel(Teuchos::VERB_EXTREME);
    Teuchos::RCP<UncoupledAggregationFactory> UCAggFact11 = Teuchos::rcp(new UncoupledAggregationFactory(/*dropFact11*/));
    UCAggFact11->SetMinNodesPerAggregate(9); // 9
    UCAggFact11->SetMaxNeighAlreadySelected(1);
    UCAggFact11->SetOrdering(MueLu::AggOptions::GRAPH);
    Teuchos::RCP<TentativePFactory> Ptent11Fact = Teuchos::rcp(new TentativePFactory(/*UCAggFact11,amalgFact11*/)); // check me
    Teuchos::RCP<TentativePFactory> P11Fact = Ptent11Fact;
    Teuchos::RCP<TransPFactory> R11Fact = Teuchos::rcp(new TransPFactory(P11Fact));
    Teuchos::RCP<NullspaceFactory> nspFact11 = Teuchos::rcp(new NullspaceFactory("Nullspace1",P11Fact));
    Teuchos::RCP<CoarseMapFactory> coarseMapFact11 = Teuchos::rcp(new CoarseMapFactory(/*UCAggFact11,nspFact11*/));
    coarseMapFact11->setStridingData(stridingInfo1);

    ///////////////////////////////////////////////////////////////////////
    // define factory manager for (1,1) block
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<FactoryManager> M11 = Teuchos::rcp(new FactoryManager());
    M11->SetFactory("A", A11Fact);
    M11->SetFactory("P", P11Fact);
    M11->SetFactory("Ptent", Ptent11Fact);
    M11->SetFactory("R", R11Fact);
    M11->SetFactory("Nullspace", nspFact11);
    M11->SetFactory("Ptent", P11Fact);
    M11->SetFactory("CoarseMap", coarseMapFact11);
    M11->SetIgnoreUserData(true);               // always use data from factories defined in factory manager

    ///////////////////////////////////////////////////////////////////////
    // create default nullspace for block 2
    ///////////////////////////////////////////////////////////////////////

    int dimNS2 = numdf;
    Teuchos::RCP<MultiVector> nspVector22 = MultiVectorFactory::Build(xA22->getRowMap(), dimNS2);

    for (int i=0; i<dimNS2; ++i) {
      Teuchos::ArrayRCP<Scalar> nsValues22 = nspVector22->getDataNonConst(i);
      int numBlocks = nsValues22.size() / dimNS2;
      for (int j=0; j< numBlocks; ++j) {
        nsValues22[j*dimNS2 + i] = 1.0;
      }
    }

    // set nullspace for block 2
    H->GetLevel(0)->Set("Nullspace2",nspVector22);

    ///////////////////////////////////////////////////////////////////////
    // set up block 2 factories
    ///////////////////////////////////////////////////////////////////////

    // use special aggregation routine which reconstructs aggregates for the
    // Lagrange multipliers using the aggregates for the displacements at the
    // contact/meshtying interface
    // keep correlation of interface nodes (displacements) and corresponding
    // Lagrange multiplier DOFs
    Teuchos::RCP<MueLu::ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > UCAggFact22 = Teuchos::rcp(new MueLu::ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(UCAggFact11, amalgFact11));

    Teuchos::RCP<AmalgamationFactory> amalgFact22 = Teuchos::rcp(new AmalgamationFactory(/*A22Fact*/));

    // use tentative prolongation operator (prolongation operator smoothing doesn't make sense since
    // the A11 block is not valid for smoothing)
    Teuchos::RCP<TentativePFactory> P22Fact = Teuchos::rcp(new TentativePFactory(/*UCAggFact22, amalgFact22*/));
    Teuchos::RCP<TransPFactory> R22Fact = Teuchos::rcp(new TransPFactory(P22Fact));
    Teuchos::RCP<NullspaceFactory> nspFact22 = Teuchos::rcp(new NullspaceFactory("Nullspace2",P22Fact));
    Teuchos::RCP<BlockedCoarseMapFactory> coarseMapFact22 = Teuchos::rcp(new BlockedCoarseMapFactory(coarseMapFact11,UCAggFact22,nspFact22));
    coarseMapFact22->setStridingData(stridingInfo2);

    ///////////////////////////////////////////////////////////////////////
    // define factory manager for (2,2) block
    ///////////////////////////////////////////////////////////////////////

    Teuchos::RCP<FactoryManager> M22 = Teuchos::rcp(new FactoryManager());
    M22->SetFactory("A", A22Fact);
    M22->SetFactory("P", P22Fact);
    M22->SetFactory("R", R22Fact);
    M22->SetFactory("Aggregates", UCAggFact22);
    M22->SetFactory("Nullspace", nspFact22);
    M22->SetFactory("Ptent", P22Fact);
    M22->SetFactory("CoarseMap", coarseMapFact22);
    M22->SetIgnoreUserData(true);               // always use data from factories defined in factory manager

    ///////////////////////////////////////////////////////////////////////
    // define block transfer operators
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<BlockedPFactory> PFact = Teuchos::rcp(new BlockedPFactory(Teuchos::null)); // use row map index base from bOp
    PFact->AddFactoryManager(M11);
    PFact->AddFactoryManager(M22);

    Teuchos::RCP<GenericRFactory> RFact = Teuchos::rcp(new GenericRFactory(PFact));

    ///////////////////////////////////////////////////////////////////////
    // define RAPFactory
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<RAPFactory> AcFact = Teuchos::rcp(new RAPFactory(PFact, RFact));
    AcFact->SetRepairZeroDiagonal(true); // repair zero diagonal entries in Ac, that are resulting from Ptent with nullspacedim > ndofspernode

    ///////////////////////////////////////////////////////////////////////
    // transfer "SlaveDofMap" to next coarser level
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cmTransFact3 = Teuchos::rcp(new MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", Ptent11Fact, MueLu::NoFactory::getRCP()));
    AcFact->AddTransferFactory(cmTransFact3);

    ///////////////////////////////////////////////////////////////////////
    // create Braess-Sarazin smoother
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<SmootherFactory> SmooFactCoarsest = GetCoarsestBraessSarazinSmootherFactory(mllist_, 0, Teuchos::null /* AFact*/);

    ///////////////////////////////////////////////////////////////////////
    // prepare factory managers
    ///////////////////////////////////////////////////////////////////////

    bool bIsLastLevel = false;
    std::vector<Teuchos::RCP<FactoryManager> > vecManager(maxLevels);
    for(int i=0; i < maxLevels; i++) {

      Teuchos::ParameterList pp(mllist_);

      // fine/intermedium level smoother
      Teuchos::RCP<SmootherFactory> SmooFactFine = GetBraessSarazinSmootherFactory(pp, i, Teuchos::null /* AFact*/);

      vecManager[i] = Teuchos::rcp(new FactoryManager());
      if(SmooFactFine != Teuchos::null)
          vecManager[i]->SetFactory("Smoother" ,  SmooFactFine);    // Hierarchy.Setup uses TOPSmootherFactory, that only needs "Smoother"
      vecManager[i]->SetFactory("CoarseSolver", SmooFactCoarsest);
      vecManager[i]->SetFactory("A", AcFact);       // same RAP factory for all levels
      vecManager[i]->SetFactory("P", PFact);        // same prolongator and restrictor factories for all levels
      vecManager[i]->SetFactory("R", RFact);        // same prolongator and restrictor factories for all levels
    }

    // use new Hierarchy::Setup routine
    if(maxLevels == 1) {
      bIsLastLevel = H->Setup(0, Teuchos::null, vecManager[0].ptr(), Teuchos::null); // 1 level "multigrid" method
    }
    else
    {
      bIsLastLevel = H->Setup(0, Teuchos::null, vecManager[0].ptr(), vecManager[1].ptr()); // first (finest) level
      for(int i=1; i < maxLevels-1; i++) { // intermedium levels
        if(bIsLastLevel == true) break;
        bIsLastLevel = H->Setup(i, vecManager[i-1].ptr(), vecManager[i].ptr(), vecManager[i+1].ptr());
      }
      if(bIsLastLevel == false) { // coarsest level
          bIsLastLevel = H->Setup(maxLevels-1, vecManager[maxLevels-2].ptr(), vecManager[maxLevels-1].ptr(), Teuchos::null);
       }
    }

    ///////////////////////////////////////////////////////////////////////
    // call setup
    ///////////////////////////////////////////////////////////////////////
    //H->Setup(M,0,maxLevels);

    // set multigrid preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > LINALG::SOLVER::MueLuContactSpPreconditioner::GetBraessSarazinSmootherFactory(const Teuchos::ParameterList & paramList, int level, const Teuchos::RCP<FactoryBase> & AFact) {
  char levelchar[11];
  sprintf(levelchar,"(level %d)",level);
  std::string levelstr(levelchar);

  if(paramList.isSublist("smoother: list " + levelstr)==false)
    return Teuchos::null;
  TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist("smoother: list " + levelstr)==false, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: no ML smoother parameter list for level. error.");

  std::string type = paramList.sublist("smoother: list " + levelstr).get<std::string>("smoother: type");
  TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: no ML smoother type for level. error.");

  const Teuchos::ParameterList smolevelsublist = paramList.sublist("smoother: list " + levelstr);

  if (type != "symmetric Gauss-Seidel") {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::ContactSPPreconditioner: Please set the ML_SMOOTHERCOARSE, ML_SMOOTHERMED and ML_SMOOTHERFINE parameters to SGS in your dat file. Other smoother options are not accepted. \n Note: In fact we're using only the ML_DAMPFINE, ML_DAMPMED, ML_DAMPCOARSE as well as the ML_SMOTIMES parameters for Braess-Sarazin.");
  }

  //smooProto = Teuchos::rcp( new MueLu::MyTrilinosSmoother<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact) );*/
  Scalar omega = smolevelsublist.get<double>("smoother: damping factor");   // Braess-Sarazin damping/scaling factor
  int sweeps   = smolevelsublist.get<int>("smoother: sweeps");

  // create SchurComp factory (SchurComplement smoother is provided by local FactoryManager)
  Teuchos::RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory(MueLu::NoFactory::getRCP(),omega));
  Teuchos::RCP<BraessSarazinSmoother> smootherPrototype = Teuchos::rcp(new BraessSarazinSmoother(sweeps,omega)); // append SC smoother information

  // define SchurComplement solver
  /*Teuchos::ParameterList SCList;
  SCList.set("relaxation: sweeps", (LO) 10);
  SCList.set("relaxation: damping factor", (SC) 0.6);
  SCList.set("relaxation: type", "Gauss-Seidel");
  Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::rcp(new TrilinosSmoother("RELAXATION",SCList,0,SFact));*/
  //Teuchos::RCP<SmootherPrototype> smoProtoSC = MueLu::GetIfpackSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("ILU", SCList,0,SFact);
  Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::rcp( new DirectSolver("Klu"/*"Umfpack"*/,Teuchos::ParameterList(),SFact) );
  Teuchos::RCP<SmootherFactory> SmooSCFact = Teuchos::rcp(new SmootherFactory(smoProtoSC));

  // setup local factory manager for SchurComplementFactory
  Teuchos::RCP<FactoryManager> MB = Teuchos::rcp(new FactoryManager());
  MB->SetFactory("A", SFact);              // SchurCompFactory as generating factory for SchurComp equation
  MB->SetFactory("Smoother", SmooSCFact);  // solver for SchurComplement equation
  MB->SetIgnoreUserData(true);
  smootherPrototype->SetFactoryManager(MB);  // add SC smoother information

  // create smoother factory
  Teuchos::RCP<SmootherFactory>   SmooFact;
  SmooFact = Teuchos::rcp( new SmootherFactory(smootherPrototype) );

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if(smolevelsublist.isParameter("smoother: pre or post")) preorpost = smolevelsublist.get<std::string>("smoother: pre or post");

  if (preorpost == "pre") {
    SmooFact->SetSmootherPrototypes(smootherPrototype, Teuchos::null);
  } else if(preorpost == "post") {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smootherPrototype);
  }

  return SmooFact;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > LINALG::SOLVER::MueLuContactSpPreconditioner::GetCoarsestBraessSarazinSmootherFactory(const Teuchos::ParameterList & paramList, int level, const Teuchos::RCP<FactoryBase> & AFact) {
  char levelchar[11];
  sprintf(levelchar,"(level %d)",level);
  std::string levelstr(levelchar);

  if(paramList.isSublist("smoother: list " + levelstr)==false)
    return Teuchos::null;
  TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist("smoother: list " + levelstr)==false, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: no ML smoother parameter list for coarsest level. error.");

  std::string type = paramList.sublist("smoother: list " + levelstr).get<std::string>("smoother: type");
  TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), MueLu::Exceptions::RuntimeError, "MueLu::ContactSpPreconditioner: no ML smoother type for level. error.");

  const Teuchos::ParameterList smolevelsublist = paramList.sublist("smoother: list " + levelstr);

  if (type != "symmetric Gauss-Seidel") {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::ContactSPPreconditioner: Please set the ML_SMOOTHERCOARSE, ML_SMOOTHERMED and ML_SMOOTHERFINE parameters to SGS in your dat file. Other smoother options are not accepted. \n Note: In fact we're using only the ML_DAMPFINE, ML_DAMPMED, ML_DAMPCOARSE as well as the ML_SMOTIMES parameters for Braess-Sarazin.");
  }

  //smooProto = Teuchos::rcp( new MueLu::MyTrilinosSmoother<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact) );*/
  Scalar omega = smolevelsublist.get<double>("smoother: damping factor");   // Braess-Sarazin damping/scaling factor
  int sweeps   = smolevelsublist.get<int>("smoother: sweeps");

  // create SchurComp factory (SchurComplement smoother is provided by local FactoryManager)
  Teuchos::RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory(MueLu::NoFactory::getRCP(),omega));
  Teuchos::RCP<BraessSarazinSmoother> smootherPrototype = Teuchos::rcp(new BraessSarazinSmoother(sweeps,omega)); // append SC smoother information

  // define SchurComplement solver
  /*Teuchos::ParameterList SCList;
  SCList.set("relaxation: sweeps", (LO) 10);
  SCList.set("relaxation: damping factor", (SC) 0.6);
  SCList.set("relaxation: type", "Gauss-Seidel");
  Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::rcp(new TrilinosSmoother("RELAXATION",SCList,0,SFact));*/
  //Teuchos::RCP<SmootherPrototype> smoProtoSC = MueLu::GetIfpackSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("ILU", SCList,0,SFact);
  Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::rcp( new DirectSolver("Klu"/*"Umfpack"*/,Teuchos::ParameterList(),SFact) );
  Teuchos::RCP<SmootherFactory> SmooSCFact = Teuchos::rcp(new SmootherFactory(smoProtoSC));

  // setup local factory manager for SchurComplementFactory
  Teuchos::RCP<FactoryManager> MB = Teuchos::rcp(new FactoryManager());
  MB->SetFactory("A", SFact);              // SchurCompFactory as generating factory for SchurComp equation
  MB->SetFactory("Smoother", SmooSCFact);  // solver for SchurComplement equation
  MB->SetIgnoreUserData(true);
  smootherPrototype->SetFactoryManager(MB);  // add SC smoother information

  // create smoother factory
  Teuchos::RCP<SmootherFactory>   SmooFact;
  SmooFact = Teuchos::rcp( new SmootherFactory(smootherPrototype) );

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if(smolevelsublist.isParameter("smoother: pre or post")) preorpost = smolevelsublist.get<std::string>("smoother: pre or post");

  if (preorpost == "pre") {
    SmooFact->SetSmootherPrototypes(smootherPrototype, Teuchos::null);
  } else if(preorpost == "post") {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smootherPrototype);
  }

  return SmooFact;
}


#endif // HAVE_MueLu

