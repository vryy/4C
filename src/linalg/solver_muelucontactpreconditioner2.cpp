/*
 * solver_muelucontactpreconditioner2.cpp
 *
 *  Created on: Aug 2, 2012
 *      Author: wiesner
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

// MueLu
#include <MueLu.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_SmootherPrototype_decl.hpp>

#include <MueLu_CoalesceDropFactory.hpp>
//#include <MueLu_UCAggregationFactory.hpp>
#include <MueLu_ExperimentalAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
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

#include <MueLu_MLParameterListInterpreter.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include <MueLu_UseShortNames.hpp>

#include <MueLu_EpetraOperator.hpp> // Aztec interface

#include "muelu_ContactAFilterFactory_decl.hpp"
#include "muelu_ContactTransferFactory_decl.hpp"
#include "muelu_ContactMapTransferFactory_decl.hpp"
#include "muelu_ContactASlaveDofFilterFactory_decl.hpp"
#include "MueLu_MyTrilinosSmoother_decl.hpp"

#include "solver_muelucontactpreconditioner2.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuContactPreconditioner2::MueLuContactPreconditioner2( FILE * outfile, Teuchos::ParameterList & mllist )
  : PreconditionerType( outfile ),
    mllist_( mllist )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuContactPreconditioner2::Setup( bool create,
                                              Epetra_Operator * matrix,
                                              Epetra_MultiVector * x,
                                              Epetra_MultiVector * b )
{
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
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Hierarchy> LINALG::SOLVER::MueLuContactPreconditioner2::SetupHierarchy(
    const Teuchos::ParameterList & params,
    const Teuchos::RCP<Operator> & A,
    const Teuchos::RCP<MultiVector> nsp)
{
//#include "MueLu_UseShortNames.hpp" // TODO don't know why this is needed here

  //std::cout << "MueLuContactPreconditoner2" << std::endl << "??????????????????????????????????????" << std::endl;

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // read in common parameters
  int maxLevels = 10;       // multigrid prameters
  int verbosityLevel = 10;  // verbosity level
  int maxCoarseSize = 50;
  int nDofsPerNode = 1;         // coalesce and drop parameters
  double agg_threshold = 0.0;   // aggregation parameters
  double agg_damping = 4/3;
  int    agg_smoothingsweeps = 1;
  int    minPerAgg = 3;       // optimal for 2d
  int    maxNbrAlreadySelected = 0;
  std::string agg_type = "Uncoupled";
  bool   bEnergyMinimization = false; // PGAMG
  if(params.isParameter("max levels")) maxLevels = params.get<int>("max levels");
  if(params.isParameter("ML output"))  verbosityLevel = params.get<int>("ML output");
  if(params.isParameter("coarse: max size")) maxCoarseSize = params.get<int>("coarse: max size");
  if(params.isParameter("PDE equations")) nDofsPerNode = params.get<int>("PDE equations");
  if(params.isParameter("aggregation: threshold"))          agg_threshold       = params.get<double>("aggregation: threshold");
  if(params.isParameter("aggregation: damping factor"))     agg_damping         = params.get<double>("aggregation: damping factor");
  if(params.isParameter("aggregation: smoothing sweeps"))   agg_smoothingsweeps = params.get<int>   ("aggregation: smoothing sweeps");
  if(params.isParameter("aggregation: type"))               agg_type            = params.get<std::string> ("aggregation: type");
  if(params.isParameter("aggregation: nodes per aggregate"))minPerAgg           = params.get<int>("aggregation: nodes per aggregate");
  if(params.isParameter("energy minimization: enable"))  bEnergyMinimization = params.get<bool>("energy minimization: enable");

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

  // build map extractor from different maps
  // note that the ordering (Master, Slave, Inner) is important to be the same overall the whole algorithm
  Teuchos::RCP<const Map> xfullmap = A->getRowMap(); // full map (MasterDofMap + SalveDofMap + InnerDofMap)
  Teuchos::RCP<Xpetra::EpetraMap> xMasterDofMap  = Teuchos::rcp(new Xpetra::EpetraMap( epMasterDofMap ));
  Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epSlaveDofMap  ));
  Teuchos::RCP<Xpetra::EpetraMap> xActiveDofMap  = Teuchos::rcp(new Xpetra::EpetraMap( epActiveDofMap ));
  //Teuchos::RCP<Xpetra::EpetraMap> xInnerDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epInnerDofMap  )); // TODO check me

  std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > > xmaps;
  xmaps.push_back(xMasterDofMap);
  xmaps.push_back(xSlaveDofMap );
  //xmaps.push_back(xInnerDofMap ); // TODO check me

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LO,GO,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LO,GO>::Build(xfullmap,xmaps);

  ///////////////////////////////////////////////////////////

  // fill hierarchy
  Teuchos::RCP<Hierarchy> hierarchy = Teuchos::rcp(new Hierarchy(A));
  hierarchy->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(eVerbLevel));
  //hierarchy->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(Teuchos::VERB_NONE));
  hierarchy->SetMaxCoarseSize(Teuchos::as<Xpetra::global_size_t>(maxCoarseSize));
  hierarchy->SetDebug(true);

  int timestep = mllist_.get<int>("time-step");
  int newtoniter = mllist_.get<int>("newton-iter");
  std::stringstream str; str << "t" << timestep << "_n" << newtoniter;
  hierarchy->SetDebugPrefix(str.str());

  ///////////////////////////////////////////////////////////

  // set fine level nullspace
  // use given fine level null space or extract pre-computed nullspace from ML parameter list
 Teuchos::RCP<MueLu::Level> Finest = hierarchy->GetLevel();  // get finest level

 Finest->Set("ActiveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xActiveDofMap));  // set map with active dofs
 Finest->Set("MasterDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xMasterDofMap));  // set map with active dofs
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

  // prepare aggCoarseStat
  
  // number of node rows
  const LocalOrdinal nDofRows = xfullmap->getNodeNumElements();
  Teuchos::ArrayRCP<MueLu::NodeState> aggStat;
  if(nDofRows > 0) aggStat = Teuchos::arcp<MueLu::NodeState>(nDofRows/nDofsPerNode);
  for(LocalOrdinal i=0; i<nDofRows; ++i) {
    aggStat[i/nDofsPerNode] = MueLu::READY;
    //if(i%333==0) aggStat[i] = MueLu::READY_1PT;
    //if((i>59 && i< 120)) aggStat[i] = MueLu::READY_1PT;
    //if((i>599 && i< 630)) aggStat[i] = MueLu::READY_1PT;
    //std::cout << "i: " << "aggStat=" << aggStat[i] << std::endl;
    GlobalOrdinal grid = xfullmap->getGlobalElement(i);
    if(xSlaveDofMap->isNodeGlobalElement(grid)) {
      aggStat[i/nDofsPerNode] = MueLu::READY_1PT;
    }
  }
  Finest->Set("coarseAggStat",aggStat);
  // create factories

  // prepare (filtered) A Factory
  //Teuchos::RCP<SingleLevelFactoryBase> segAFact = Teuchos::rcp(new MueLu::ContactAFilterFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("A", NULL, map_extractor));

  // for the Jacobi/SGS smoother we wanna change the input matrix A and set Dirichlet bcs for the (active?) slave dofs
  Teuchos::RCP<FactoryBase> slaveDcAFact = Teuchos::rcp(new MueLu::ContactASlaveDofFilterFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>());
  Finest->Keep("A",slaveDcAFact.get()); // TODO remove me? no, since it is used for the level smoothers!!!

  // Coalesce and drop factory with constant number of Dofs per freedom
  Teuchos::RCP<CoalesceDropFactory> dropFact = Teuchos::rcp(new CoalesceDropFactory(/*segAFact*/Teuchos::null,Teuchos::null));

  // aggregation factory
  //Teuchos::RCP<UCAggregationFactory> UCAggFact = Teuchos::rcp(new UCAggregationFactory(dropFact));
  Teuchos::RCP<ExperimentalAggregationFactory> UCAggFact = Teuchos::rcp(new ExperimentalAggregationFactory(dropFact));
  UCAggFact->SetMinNodesPerAggregate(minPerAgg); //TODO should increase if run anything other than 1D
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);

  // transfer operators (PG-AMG)
  Teuchos::RCP<PFactory> PFact = Teuchos::rcp(new TentativePFactory(UCAggFact,Teuchos::null,Teuchos::null,Teuchos::null));
  Teuchos::RCP<RFactory> RFact  = Teuchos::rcp( new TransPFactory(PFact) );
  /*Teuchos::RCP<PFactory> PtentFact = Teuchos::rcp(new TentativePFactory(UCAggFact,Teuchos::null,Teuchos::null,Teuchos::null));
  Teuchos::RCP<PFactory> PFact  = Teuchos::rcp( new PgPFactory(PtentFact) );
  Teuchos::RCP<RFactory> RFact  = Teuchos::rcp( new GenericRFactory(PFact) );*/

  // define nullspace factory AFTER tentative PFactory (that generates the nullspace for the coarser levels)
  // use same nullspace factory for all multigrid levels
  // therefor we have to create one instance of NullspaceFactory and use it
  // for all FactoryManager objects (note: here, we have one FactoryManager object per level)
  Teuchos::RCP<NullspaceFactory> nspFact = Teuchos::rcp(new NullspaceFactory("Nullspace",PFact /* TODO tentativePFactory */));
  nspFact->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(eVerbLevel));
  //Teuchos::RCP<NullspaceFactory> nspFact = Teuchos::rcp(new NullspaceFactory("Nullspace",PtentFact /* TODO tentativePFactory */));

  // RAP factory with inter-level transfer of segregation block information (map extractor)
  Teuchos::RCP<RAPFactory> AcFact = Teuchos:: rcp( new RAPFactory(PFact, RFact) );
  AcFact->setVerbLevel(Teuchos::VERB_HIGH);
  AcFact->SetRepairZeroDiagonal(true); // repair zero diagonal entries in Ac

  // write out aggregates
  Teuchos::RCP<MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > aggExpFact = Teuchos::rcp(new MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>("aggs_level%LEVELID_proc%PROCID.out",UCAggFact.get(), dropFact.get(),NULL/*amalgFact*/));
  AcFact->AddTransferFactory(aggExpFact);

  //Teuchos::RCP<MueLu::ContactTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cTransFact = Teuchos::rcp(new MueLu::ContactTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(PFact));
  //AcFact->AddTransferFactory(cTransFact);

  // transfer maps to coarser grids

  Teuchos::RCP<MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cmTransFact = Teuchos::rcp(new MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("ActiveDofMap", PFact, MueLu::NoFactory::getRCP()));
  AcFact->AddTransferFactory(cmTransFact);
  Teuchos::RCP<MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cmTransFact2 = Teuchos::rcp(new MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("MasterDofMap", PFact, MueLu::NoFactory::getRCP()));
  AcFact->AddTransferFactory(cmTransFact2);
  Teuchos::RCP<MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cmTransFact3 = Teuchos::rcp(new MueLu::ContactMapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", PFact, MueLu::NoFactory::getRCP()));
  AcFact->AddTransferFactory(cmTransFact3);

  Teuchos::RCP<MueLu::AggStatTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> > aggStatFact = Teuchos::rcp(new MueLu::AggStatTransferFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps>("coarseAggStat",UCAggFact));
  AcFact->AddTransferFactory(aggStatFact);


  // setup smoothers
  Teuchos::RCP<SmootherFactory> coarsestSmooFact;
  coarsestSmooFact = MueLu::MLParameterListInterpreter<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetCoarsestSolverFactory(params);

  ///////////////////////////////////////////////////




  Finest->Set("A",A);
  Finest->Keep("Aggregates",UCAggFact.get());

  ////////////////////////////////////

  // prepare factory managers

  bool bIsLastLevel = false;
  std::vector<Teuchos::RCP<FactoryManager> > vecManager(maxLevels);
  for(int i=0; i < maxLevels; i++) {
    //params.set("smoother: pre or post","pre"); // only pre-smoothing
    Teuchos::ParameterList pp(params);
    //pp.set("smoother: pre or post","pre");

    //Teuchos::RCP<SmootherFactory> SmooFactFine = MueLu::MLParameterListInterpreter<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetSmootherFactory(pp, i, slaveDcAFact);

    //////////////////////////////////////
    Teuchos::RCP<SmootherPrototype> smooProto;
    std::string ifpackType;
    Teuchos::ParameterList ifpackList;
    Teuchos::RCP<SmootherFactory> SmooFactFine;

    ifpackType = "RELAXATION";
    ifpackList.set<int>("relaxation: sweeps", 1);
    ifpackList.set("relaxation: damping factor", 0.6);
    ifpackList.set("relaxation: type", "Jacobi");

    smooProto = Teuchos::rcp( new MueLu::MyTrilinosSmoother<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, slaveDcAFact) );
    SmooFactFine = Teuchos::rcp( new SmootherFactory(smooProto) );

    //////////////////////////////////////

    vecManager[i] = Teuchos::rcp(new FactoryManager());
    if(SmooFactFine != Teuchos::null)
        vecManager[i]->SetFactory("Smoother" ,  SmooFactFine);    // Hierarchy.Setup uses TOPSmootherFactory, that only needs "Smoother"
    vecManager[i]->SetFactory("CoarseSolver", coarsestSmooFact);
    vecManager[i]->SetFactory("Aggregates", UCAggFact);
    vecManager[i]->SetFactory("Graph", dropFact);
    vecManager[i]->SetFactory("DofsPerNode", dropFact);
    vecManager[i]->SetFactory("A", AcFact);       // same RAP factory
    vecManager[i]->SetFactory("P", PFact);    // same prolongator and restrictor factories
    vecManager[i]->SetFactory("Ptent", PFact);    // same prolongator and restrictor factories
    vecManager[i]->SetFactory("R", RFact);    // same prolongator and restrictor factories
    vecManager[i]->SetFactory("Nullspace", nspFact); // use same nullspace factory throughout all multigrid levels
  }

  // use new Hierarchy::Setup routine
  if(maxLevels == 1) {
          bIsLastLevel = hierarchy->Setup(0, Teuchos::null, vecManager[0].ptr(), Teuchos::null);
  }
  else
  {
          bIsLastLevel = hierarchy->Setup(0, Teuchos::null, vecManager[0].ptr(), vecManager[1].ptr()); // true, false because first level
          for(int i=1; i < maxLevels-1; i++) {
                if(bIsLastLevel == true) break;
                bIsLastLevel = hierarchy->Setup(i, vecManager[i-1].ptr(), vecManager[i].ptr(), vecManager[i+1].ptr());
          }
          if(bIsLastLevel == false) {
                bIsLastLevel = hierarchy->Setup(maxLevels-1, vecManager[maxLevels-2].ptr(), vecManager[maxLevels-1].ptr(), Teuchos::null);
          }
  }


  return hierarchy;
}



#endif

