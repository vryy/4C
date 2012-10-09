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


    // fix null space for "Inverse1"
    //      {
    //        const Epetra_Map& oldmap = A->FullRowMap();
    //        const Epetra_Map& newmap = A->Matrix(0,0).EpetraMatrix()->RowMap();
    //        LINALG::Solver::FixMLNullspace("Inverse1",oldmap, newmap, params_.sublist("Inverse1"));
    //      }

    // adapt null space for constraint equations
    /*Teuchos::ParameterList& inv2 = mllist_.sublist("Inverse2");
    if(inv2.isSublist("ML Parameters"))
    {
      // Schur complement system (1 degree per "node") -> standard nullspace
      inv2.sublist("ML Parameters").set("PDE equations",1);
      inv2.sublist("ML Parameters").set("null space: dimension",1);
      const int plength = (*A)(1,1).RowMap().NumMyElements();
      Teuchos::RCP<std::vector<double> > pnewns = Teuchos::rcp(new std::vector<double>(plength,1.0));
      //TODO: vector<double> has zero length for particular cases (e.g. no Lagrange multiplier on this processor)
      //      -> RCP for the null space is set to NULL in Fedora 12 -> dserror
      //      -> RCP points to a random memory field in Fedora 8 -> RCP for null space is not NULL
      // Temporary work around (ehrl, 21.12.11):
      // In the case of plength=0 the vector<double> is rescaled (size 0 -> size 1, initial value 0) in order to avoid problems with ML
      // (ML expects an RCP for the null space != NULL)
      if (plength==0)
        pnewns->resize(1,0.0);
      inv2.sublist("ML Parameters").set("null space: vectors",&((*pnewns)[0]));
      inv2.sublist("ML Parameters").remove("nullspace",false);
      inv2.sublist("Michael's secret vault").set<Teuchos::RCP<std::vector<double> > >("pressure nullspace",pnewns);
    }*/

     // create a Teuchos::Comm from EpetraComm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::toXpetra(A->RangeMap(0).Comm());

    // TODO build strided maps (using the maps from block matrix A...

    Teuchos::RCP<const Map> fullrangemap = Teuchos::rcp(new Xpetra::EpetraMap(Teuchos::rcpFromRef(A->FullRangeMap())));

    Teuchos::RCP<CrsMatrix> xA11 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0,0).EpetraMatrix()));
    Teuchos::RCP<CrsMatrix> xA12 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0,1).EpetraMatrix()));
    Teuchos::RCP<CrsMatrix> xA21 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1,0).EpetraMatrix()));
    Teuchos::RCP<CrsMatrix> xA22 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1,1).EpetraMatrix()));

    ///////////////////// EXPERIMENTAL
    std::vector<size_t> stridingInfo1;
    stridingInfo1.push_back(numdf);
    Teuchos::RCP<Xpetra::StridedEpetraMap> strMap1 = Teuchos::rcp(new Xpetra::StridedEpetraMap(Teuchos::rcpFromRef(A->Matrix(0,0).EpetraMatrix()->RowMap()), stridingInfo1, -1 /* stridedBlock */, 0 /*globalOffset*/));
    std::vector<size_t> stridingInfo2;
    stridingInfo2.push_back(1);
    Teuchos::RCP<Xpetra::StridedEpetraMap> strMap2 = Teuchos::rcp(new Xpetra::StridedEpetraMap(Teuchos::rcpFromRef(A->Matrix(1,1).EpetraMatrix()->RowMap()), stridingInfo2, -1 /* stridedBlock */, 0 /*globalOffset*/));
    std::cout << *strMap1 << std::endl;
    std::cout << *strMap2 << std::endl;
    ///////////////////// EXPERIMENTAL

    std::vector<Teuchos::RCP<const Map> > xmaps;
    //xmaps.push_back(xA11->getRowMap()); // TODO introduce strided maps
    //xmaps.push_back(xA22->getRowMap());
    xmaps.push_back(strMap1);
    xmaps.push_back(strMap2);

    std::cout << "length of A11 block: " << xA11->getRowMap()->getGlobalNumElements() << std::endl;

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



    std::cout << "PDE equations: " << numdf << std::endl;
    std::cout << "nullspace dims: " << dimns << std::endl;

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
    }

    // create Hierarchy
    Teuchos::RCP<Hierarchy> H = rcp(new Hierarchy());
    H->setDefaultVerbLevel(Teuchos::VERB_EXTREME);
    H->SetMaxCoarseSize(100); // TODO fix me
    H->GetLevel(0)->Set("A",Teuchos::rcp_dynamic_cast<Matrix>(bOp));
    H->GetLevel(0)->Set("Nullspace1",nspVector11);
    H->GetLevel(0)->Set("coarseAggStat",aggStat);

    Teuchos::RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));
    Teuchos::RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 1));

    // set up block 11
    Teuchos::RCP<AmalgamationFactory> amalgFact11 = Teuchos::rcp(new AmalgamationFactory(A11Fact));
    amalgFact11->setDefaultVerbLevel(Teuchos::VERB_EXTREME);
    Teuchos::RCP<CoalesceDropFactory> dropFact11 = Teuchos::rcp(new CoalesceDropFactory(A11Fact,amalgFact11));
    dropFact11->setDefaultVerbLevel(Teuchos::VERB_EXTREME);
    Teuchos::RCP<UncoupledAggregationFactory> UCAggFact11 = Teuchos::rcp(new UncoupledAggregationFactory(dropFact11));
    UCAggFact11->SetMinNodesPerAggregate(9);
    UCAggFact11->SetMaxNeighAlreadySelected(1);
    UCAggFact11->SetOrdering(MueLu::AggOptions::NATURAL);
    Teuchos::RCP<TentativePFactory> P11Fact = Teuchos::rcp(new TentativePFactory(UCAggFact11,amalgFact11)); // check me
    P11Fact->setStridingData(stridingInfo1);
    //P11Fact->setStridedBlockId(0); // declare this P11Fact to be the transfer operator for the velocity dofs
    Teuchos::RCP<TransPFactory> R11Fact = Teuchos::rcp(new TransPFactory(P11Fact));
    Teuchos::RCP<NullspaceFactory> nspFact11 = Teuchos::rcp(new NullspaceFactory("Nullspace1",P11Fact));

    //////////////////////////////// define factory manager for (1,1) block
    Teuchos::RCP<FactoryManager> M11 = Teuchos::rcp(new FactoryManager());
    M11->SetFactory("A", A11Fact);
    M11->SetFactory("P", P11Fact);
    M11->SetFactory("R", R11Fact);
    M11->SetFactory("Nullspace", nspFact11);
    M11->SetFactory("Ptent", P11Fact);
    M11->SetIgnoreUserData(true);               // always use data from factories defined in factory manager


    Teuchos::RCP<MultiVector> nspVector22 = MultiVectorFactory::Build(xA22->getRowMap(), 1);  // this is a 2D standard null space
    Teuchos::ArrayRCP<Scalar> nsValues22 = nspVector22->getDataNonConst(0);
    for (int j=0; j< nsValues22.size(); ++j) {
      nsValues22[j] = 1.0;
    }

    H->GetLevel(0)->Set("Nullspace2",nspVector22);

    // use TentativePFactory
    Teuchos::RCP<AmalgamationFactory> amalgFact22 = Teuchos::rcp(new AmalgamationFactory(A22Fact));
    Teuchos::RCP<TentativePFactory> P22Fact = Teuchos::rcp(new TentativePFactory(UCAggFact11, amalgFact22));
    P22Fact->setStridingData(stridingInfo2);
    //P22Fact->setStridedBlockId(0); // declare this P22Fact to be the transfer operator for the pressure dofs

    Teuchos::RCP<TransPFactory> R22Fact = Teuchos::rcp(new TransPFactory(P22Fact));

    Teuchos::RCP<NullspaceFactory> nspFact22 = Teuchos::rcp(new NullspaceFactory("Nullspace2",P22Fact));

    //////////////////////////////// define factory manager for (2,2) block
    Teuchos::RCP<FactoryManager> M22 = Teuchos::rcp(new FactoryManager());
    M22->SetFactory("A", A22Fact);
    M22->SetFactory("P", P22Fact);
    M22->SetFactory("R", R22Fact);
    M22->SetFactory("Aggregates", UCAggFact11);
    M22->SetFactory("Nullspace", nspFact22);
    M22->SetFactory("Ptent", P22Fact);
    M22->SetIgnoreUserData(true);               // always use data from factories defined in factory manager

    /////////////////////////////////////////// define blocked transfer ops
    Teuchos::RCP<BlockedPFactory> PFact = Teuchos::rcp(new BlockedPFactory(Teuchos::null)); // use row map index base from bOp
    PFact->AddFactoryManager(M11);
    PFact->AddFactoryManager(M22);

    Teuchos::RCP<GenericRFactory> RFact = Teuchos::rcp(new GenericRFactory(PFact));

    Teuchos::RCP<RAPFactory> AcFact = Teuchos::rcp(new RAPFactory(PFact, RFact));

    // create Braess-Sarazin smoother
    Scalar omega = 1.3;
    Teuchos::RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory(MueLu::NoFactory::getRCP(),omega));
    Teuchos::RCP<BraessSarazinSmoother> smootherPrototype = Teuchos::rcp(new BraessSarazinSmoother(3,omega)); // append SC smoother information

    Teuchos::RCP<SmootherFactory> smootherFact = Teuchos::rcp(new SmootherFactory(smootherPrototype));

    // SchurComplement smoother
    Teuchos::ParameterList SCList;
    SCList.set("relaxation: sweeps", (LO) 3);
    SCList.set("relaxation: damping factor", (SC) 1.0);
    SCList.set("relaxation: type", "Gauss-Seidel");
    Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::rcp(new TrilinosSmoother("RELAXATION",SCList,0,SFact));
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

  dserror("EXIT");
#if 0
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>( matrix );
    if ( A==NULL )
      dserror( "CrsMatrix expected" );

    // free old matrix first
    P_       = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

    // see whether we use standard ml or our own mlapi operator
    //const bool domuelupreconditioner = mllist_.get<bool>("LINALG::MueLu_Preconditioner",false);

    // wrap Epetra_CrsMatrix to Xpetra::Operator for use in MueLu
    Teuchos::RCP<Xpetra::CrsMatrix<SC,LO,GO,NO,LMO > > mueluA  = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
    Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO,LMO> >   mueluOp = Teuchos::rcp(new Xpetra::CrsOperator<SC,LO,GO,NO,LMO>(mueluA));

    // prepare nullspace vector for MueLu
    int numdf = mllist_.get<int>("PDE equations",-1);
    int dimns = mllist_.get<int>("null space: dimension",-1);
    if(dimns == -1 || numdf == -1) dserror("Error: PDE equations or null space dimension wrong.");
    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap = mueluA->getRowMap();

    Teuchos::RCP<MultiVector> nspVector = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
    Teuchos::RCP<std::vector<double> > nsdata = mllist_.get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

    for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
        Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
        const size_t myLength = nspVector->getLocalLength();
        for(size_t j=0; j<myLength; j++) {
                nspVectori[j] = (*nsdata)[i*myLength+j];
        }
    }

    // remove unsupported flags
    mllist_.remove("aggregation: threshold",false); // no support for aggregation: threshold TODO

    // Setup MueLu Hierarchy
    //Teuchos::RCP<Hierarchy> H = MLInterpreter::Setup(mllist_, mueluOp, nspVector);
    Teuchos::RCP<Hierarchy> H = SetupHierarchy(mllist_, mueluOp, nspVector);

    // set preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

  }
#endif
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Hierarchy> LINALG::SOLVER::MueLuContactSpPreconditioner::SetupHierarchy(
    const Teuchos::ParameterList & params,
    const Teuchos::RCP<Matrix> & A,
    const Teuchos::RCP<MultiVector> nsp)
{

  //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // read in common parameters
  int maxLevels = 10;       // multigrid prameters
  int verbosityLevel = 10;  // verbosity level
  int maxCoarseSize = 50;
  int nDofsPerNode = 1;         // coalesce and drop parameters
  //double agg_threshold = 0.0;   // aggregation parameters
  double agg_damping = 4/3;
  //int    agg_smoothingsweeps = 1;
  int    minPerAgg = 3;       // optimal for 2d
  int    maxNbrAlreadySelected = 0;
  std::string agg_type = "Uncoupled";
  //bool   bEnergyMinimization = false; // PGAMG
  if(params.isParameter("max levels")) maxLevels = params.get<int>("max levels");
  if(params.isParameter("ML output"))  verbosityLevel = params.get<int>("ML output");
  if(params.isParameter("coarse: max size")) maxCoarseSize = params.get<int>("coarse: max size");
  if(params.isParameter("PDE equations")) nDofsPerNode = params.get<int>("PDE equations");
  //if(params.isParameter("aggregation: threshold"))          agg_threshold       = params.get<double>("aggregation: threshold");
  if(params.isParameter("aggregation: damping factor"))     agg_damping         = params.get<double>("aggregation: damping factor");
  //if(params.isParameter("aggregation: smoothing sweeps"))   agg_smoothingsweeps = params.get<int>   ("aggregation: smoothing sweeps");
  if(params.isParameter("aggregation: type"))               agg_type            = params.get<std::string> ("aggregation: type");
  if(params.isParameter("aggregation: nodes per aggregate"))minPerAgg           = params.get<int>("aggregation: nodes per aggregate");
  //if(params.isParameter("energy minimization: enable"))  bEnergyMinimization = params.get<bool>("energy minimization: enable");
#if 0
  // set DofsPerNode in A operator
  A->SetFixedBlockSize(nDofsPerNode);

  // translate verbosity parameter
  Teuchos::EVerbosityLevel eVerbLevel = Teuchos::VERB_NONE;
  if(verbosityLevel == 0)  eVerbLevel = Teuchos::VERB_NONE;
  if(verbosityLevel > 0 )  eVerbLevel = Teuchos::VERB_LOW;
  if(verbosityLevel > 4 )  eVerbLevel = Teuchos::VERB_MEDIUM;
  if(verbosityLevel > 7 )  eVerbLevel = Teuchos::VERB_HIGH;
  if(verbosityLevel > 9 )  eVerbLevel = Teuchos::VERB_EXTREME;

  // extract additional maps from parameter list
  // these maps are provided by the STR::TimInt::PrepareContactMeshtying routine, that
  // has access to the contact manager class
  Teuchos::RCP<Epetra_Map> epMasterDofMap = params.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::MasterDofMap");
  Teuchos::RCP<Epetra_Map> epSlaveDofMap  = params.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::SlaveDofMap");
  Teuchos::RCP<Epetra_Map> epActiveDofMap = params.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::ActiveDofMap");
  //Teuchos::RCP<Epetra_Map> epInnerDofMap  = params.get<Teuchos::RCP<Epetra_Map> >("LINALG::SOLVER::MueLu_ContactPreconditioner::InnerDofMap"); // TODO check me

  //std::cout << *epActiveDofMap << std::endl;

  // build map extractor from different maps
  // note that the ordering (Master, Slave, Inner) is important to be the same overall the whole algorithm
  Teuchos::RCP<const Map> xfullmap = A->getRowMap(); // full map (MasterDofMap + SalveDofMap + InnerDofMap)
  //Teuchos::RCP<Xpetra::EpetraMap> xMasterDofMap  = Teuchos::rcp(new Xpetra::EpetraMap( epMasterDofMap ));
  Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epSlaveDofMap  ));
  Teuchos::RCP<Xpetra::EpetraMap> xActiveDofMap  = Teuchos::rcp(new Xpetra::EpetraMap( epActiveDofMap ));
  //Teuchos::RCP<Xpetra::EpetraMap> xInnerDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epInnerDofMap  )); // TODO check me

  /*std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > > xmaps;
  xmaps.push_back(xMasterDofMap);
  xmaps.push_back(xSlaveDofMap );
  //xmaps.push_back(xInnerDofMap ); // TODO check me

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xfullmap,xmaps);*/

  ///////////////////////////////////////////////////////////

  // fill hierarchy
  Teuchos::RCP<Hierarchy> hierarchy = Teuchos::rcp(new Hierarchy(A));
  hierarchy->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(eVerbLevel));
  hierarchy->SetMaxCoarseSize(Teuchos::as<Xpetra::global_size_t>(maxCoarseSize));
  //hierarchy->SetDebug(true);

  /*int timestep = mllist_.get<int>("time-step");
  int newtoniter = mllist_.get<int>("newton-iter");
  std::stringstream str; str << "t" << timestep << "_n" << newtoniter;
  hierarchy->SetDebugPrefix(str.str());*/

  ///////////////////////////////////////////////////////////

  // set fine level nullspace
  // use given fine level null space or extract pre-computed nullspace from ML parameter list
 Teuchos::RCP<MueLu::Level> Finest = hierarchy->GetLevel();  // get finest level

 Finest->Set("A",A);

 Finest->Set("ActiveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xActiveDofMap));  // set map with active dofs
 //Finest->Set("MasterDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xMasterDofMap));  // set map with active dofs
 Finest->Set("SlaveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xSlaveDofMap));  // set map with active dofs

  if (nsp != Teuchos::null) {
    Finest->Set("Nullspace",nsp);                       // set user given null space
  } else {
    std::string type = "";
    if(params.isParameter("null space: type")) type = params.get<std::string>("null space: type");
    if(type != "pre-computed") dserror("MueLu::Interpreter: no valid nullspace (no pre-computed null space). error.");
    int dimns = -1;
    if(params.isParameter("null space: dimension")) dimns = params.get<int>("null space: dimension");
    if(dimns == -1) dserror( "MueLu::Interpreter: no valid nullspace (nullspace dim = -1). error.");

    const Teuchos::RCP<const Map> rowMap = A->getRowMap();
    Teuchos::RCP<MultiVector> nspVector = MultiVectorFactory::Build(rowMap,dimns,true);
    double* nsdata = NULL;
    if(params.isParameter("null space: vectors")) nsdata = params.get<double*>("null space: vectors");
    if(nsdata == NULL) dserror("MueLu::Interpreter: no valid nullspace (nsdata = NULL). error.");

    for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
      Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
      const size_t myLength = nspVector->getLocalLength();
      for(size_t j=0; j<myLength; j++) {
        nspVectori[j] = nsdata[i*myLength+j];
      }

    }
    Finest->Set("Nullspace",nspVector);                       // set user given null space
  }

  ///////////////////////////////////////////////////////////////////////
  // special aggregation strategy
  //   - use 1pt aggregates for slave nodes
  ///////////////////////////////////////////////////////////////////////

  // number of node rows
  const LocalOrdinal nDofRows = xfullmap->getNodeNumElements();

  // prepare aggCoarseStat
  // TODO rebuild node-based map
  // still problematich for reparitioning
  Teuchos::ArrayRCP<unsigned int> aggStat;
  if(nDofRows > 0) aggStat = Teuchos::arcp<unsigned int>(nDofRows/nDofsPerNode);
  for(LocalOrdinal i=0; i<nDofRows; ++i) {
    aggStat[i/nDofsPerNode] = 0;
    GlobalOrdinal grid = xfullmap->getGlobalElement(i);
    if(xActiveDofMap->isNodeGlobalElement(grid)) {
      aggStat[i/nDofsPerNode] |= MueLu::NODEONEPT;
      //std::cout << "node " << i/nDofsPerNode << " is marked as 1PT node" << std::endl;
    }
    else if(xSlaveDofMap->isNodeGlobalElement(grid)) {
      aggStat[i/nDofsPerNode] |= MueLu::NODEONEPT; // MueLu::NODESMALLAGG;
      //std::cout << "node " << i/nDofsPerNode << " is marked as SMALL AGG node" << std::endl;
    }
  }

  /*for(LocalOrdinal i=0; i<nDofRows; ++i) {
    std::cout << i << ": " << aggStat[i/nDofsPerNode] << std::endl;
  }*/

  Finest->Set("coarseAggStat",aggStat);

  ///////////////////////////////////////////////////////////////////////
  // Segregation Factory for building aggregates that do not overlap
  // contact interface
  //   - currently not used, but would make sense in future
  ///////////////////////////////////////////////////////////////////////

  // prepare (filtered) A Factory
  //Teuchos::RCP<SingleLevelFactoryBase> segAFact = Teuchos::rcp(new MueLu::ContactAFilterFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("A", NULL, map_extractor));

  ///////////////////////////////////////////////////////////////////////
  // ContactASlaveDofFilterFactory
  //   create a matrix A with artificial Dirichlet Bcs conditions on slave
  //   dofs -> avoid zeros on diagonal
  //   This is needed for level smoothers as well as prolongator smoothing
  ///////////////////////////////////////////////////////////////////////

  // for the Jacobi/SGS smoother we wanna change the input matrix A and set Dirichlet bcs for the (active?) slave dofs
  Teuchos::RCP<FactoryBase> slaveDcAFact = Teuchos::rcp(new MueLu::ContactASlaveDofFilterFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>());
  Finest->Keep("A",slaveDcAFact.get()); // do not forget to keep A for level smoothers!

  // Coalesce and drop factory with constant number of Dofs per freedom
  // note: coalescing based on original matrix A
  Teuchos::RCP<CoalesceDropFactory> dropFact = Teuchos::rcp(new CoalesceDropFactory());

  // aggregation factory
  //Teuchos::RCP<UCAggregationFactory> UCAggFact = Teuchos::rcp(new UCAggregationFactory(dropFact));
  Teuchos::RCP<ExperimentalAggregationFactory> UCAggFact = Teuchos::rcp(new ExperimentalAggregationFactory(dropFact));
  UCAggFact->SetMinNodesPerAggregate(minPerAgg);
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  //UCAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);

  Teuchos::RCP<PFactory> PFact;
  Teuchos::RCP<RFactory> RFact;

  Teuchos::RCP<PFactory> PtentFact = Teuchos::rcp(new TentativePFactory(UCAggFact));

  // choose either nonsmoothed transfer operators or
  // PG-AMG smoothed aggregation transfer operators
  // note:
  //  - SA-AMG is not working properly (probably due to problematic Dinv scaling with zeros on diagonal) TODO handling of zeros on diagonal in SaPFactory
  //  - PG-AMG has some special handling for zeros on diagonal (avoid NaNs)
  //    avoid local damping factory omega==1 -> oversmoothing, leads to zero rows in P
  //    use matrix A with artificial Dirichlet bcs for prolongator smoothing
  // if agg_damping == 0.0 -> PA-AMG else PG-AMG
  if (agg_damping == 0.0) {
    // tentative prolongation operator (PA-AMG)
    PFact = PtentFact;
    RFact = Teuchos::rcp( new TransPFactory(PFact) );
  } else {
    // Petrov Galerkin PG-AMG smoothed aggregation (energy minimization in ML)
    PFact  = Teuchos::rcp( new PgPFactory(PtentFact,slaveDcAFact) ); // use slaveDcAFact for prolongator smoothing
    RFact  = Teuchos::rcp( new GenericRFactory() );
  }

  // define nullspace factory AFTER tentative PFactory (that generates the nullspace for the coarser levels)
  // use same nullspace factory for all multigrid levels
  // therefor we have to create one instance of NullspaceFactory and use it
  // for all FactoryManager objects (note: here, we have one FactoryManager object per level)
  Teuchos::RCP<NullspaceFactory> nspFact = Teuchos::rcp(new NullspaceFactory("Nullspace",PtentFact));

  // RAP factory with inter-level transfer of segregation block information (map extractor)
  Teuchos::RCP<RAPFactory> AcFact = Teuchos:: rcp( new RAPFactory(PFact, RFact) );
  //AcFact->setVerbLevel(Teuchos::VERB_HIGH);
  AcFact->SetRepairZeroDiagonal(true); // repair zero diagonal entries in Ac, that are resulting from Ptent with nullspacedim > ndofspernode


  //Params().sublist("MueLu (Contact2) Parameters").set("time-step",Params().get<int>("time-step"));
  //Params().sublist("MueLu (Contact2) Parameters").set("newton-iter",Params().get<int>("newton-iter"));

  // write out aggregates
  /*int timestep = mllist_.get<int>("time-step");
  int newtoniter = mllist_.get<int>("newton-iter");
  std::stringstream str;
  str << "aggs";
  str << "t" << timestep << "_n" << newtoniter;
  str << "_level%LEVELID_proc%PROCID.out";*/

  Teuchos::RCP<MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > aggExpFact = Teuchos::rcp(new MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("aggs_level%LEVELID_proc%PROCID.out"/*str.str()*/,UCAggFact.get(), dropFact.get(),NULL/*amalgFact*/));
  AcFact->AddTransferFactory(aggExpFact);

  // transfer maps to coarser grids
  Teuchos::RCP<MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cmTransFact = Teuchos::rcp(new MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("ActiveDofMap", PtentFact, MueLu::NoFactory::getRCP()));
  AcFact->AddTransferFactory(cmTransFact);
  /*Teuchos::RCP<MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cmTransFact2 = Teuchos::rcp(new MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("MasterDofMap", PtentFact, MueLu::NoFactory::getRCP()));
  AcFact->AddTransferFactory(cmTransFact2);*/
  Teuchos::RCP<MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cmTransFact3 = Teuchos::rcp(new MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", PtentFact, MueLu::NoFactory::getRCP()));
  AcFact->AddTransferFactory(cmTransFact3);

  // transfer aggregate status to next coarser level (-> special aggregation strategy)
  Teuchos::RCP<MueLu::AggStatTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> > aggStatFact = Teuchos::rcp(new MueLu::AggStatTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps>("coarseAggStat",UCAggFact));
  AcFact->AddTransferFactory(aggStatFact);

  ///////////////////////////////////////////////////////////////////////
  // setup coarse level smoothers/solvers
  ///////////////////////////////////////////////////////////////////////

  // coarse level smoother/solver
  Teuchos::RCP<SmootherFactory> coarsestSmooFact;
  coarsestSmooFact = MueLu::MLParameterListInterpreter<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetCoarsestSolverFactory(params);

  ///////////////////////////////////////////////////////////////////////
  // prepare factory managers
  ///////////////////////////////////////////////////////////////////////

  bool bIsLastLevel = false;
  std::vector<Teuchos::RCP<FactoryManager> > vecManager(maxLevels);
  for(int i=0; i < maxLevels; i++) {
    //params.set("smoother: pre or post","pre"); // only pre-smoothing
    Teuchos::ParameterList pp(params);
    //pp.set("smoother: pre or post","pre");

    // fine/intermedium level smoother
    Teuchos::RCP<SmootherFactory> SmooFactFine = GetContactSmootherFactory(pp, i, slaveDcAFact);

    vecManager[i] = Teuchos::rcp(new FactoryManager());
    if(SmooFactFine != Teuchos::null)
        vecManager[i]->SetFactory("Smoother" ,  SmooFactFine);    // Hierarchy.Setup uses TOPSmootherFactory, that only needs "Smoother"
    vecManager[i]->SetFactory("CoarseSolver", coarsestSmooFact);
    vecManager[i]->SetFactory("Aggregates", UCAggFact);
    vecManager[i]->SetFactory("Graph", dropFact);
    vecManager[i]->SetFactory("DofsPerNode", dropFact);
    vecManager[i]->SetFactory("A", AcFact);       // same RAP factory for all levels
    vecManager[i]->SetFactory("P", PFact);        // same prolongator and restrictor factories for all levels
    vecManager[i]->SetFactory("Ptent", PtentFact);// same prolongator and restrictor factories for all levels
    vecManager[i]->SetFactory("R", RFact);        // same prolongator and restrictor factories for all levels
    vecManager[i]->SetFactory("Nullspace", nspFact); // use same nullspace factory throughout all multigrid levels
  }

  // use new Hierarchy::Setup routine
  if(maxLevels == 1) {
    bIsLastLevel = hierarchy->Setup(0, Teuchos::null, vecManager[0].ptr(), Teuchos::null); // 1 level "multigrid" method
  }
  else
  {
    bIsLastLevel = hierarchy->Setup(0, Teuchos::null, vecManager[0].ptr(), vecManager[1].ptr()); // first (finest) level
    for(int i=1; i < maxLevels-1; i++) { // intermedium levels
      if(bIsLastLevel == true) break;
      bIsLastLevel = hierarchy->Setup(i, vecManager[i-1].ptr(), vecManager[i].ptr(), vecManager[i+1].ptr());
    }
    if(bIsLastLevel == false) { // coarsest level
        bIsLastLevel = hierarchy->Setup(maxLevels-1, vecManager[maxLevels-2].ptr(), vecManager[maxLevels-1].ptr(), Teuchos::null);
     }
  }

  return hierarchy;
#else
  return Teuchos::null;
#endif

}


#endif // HAVE_MueLu

