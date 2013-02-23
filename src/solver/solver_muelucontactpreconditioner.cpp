/*
 * solver_muelucontactpreconditioner.cpp
 *
 *  Created on: Aug 2, 2012
 *      Author: wiesner
 */


#ifdef HAVE_MueLu
#ifdef HAVE_EXPERIMENTAL_MueLu

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
#include <MueLu_DirectSolver.hpp>
#include <MueLu_SmootherPrototype_decl.hpp>

#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_PgPFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_NullspaceFactory.hpp>
//#include <MueLu_SegregationATransferFactory.hpp> // TODO remove me
#include <MueLu_Aggregates.hpp>
#include <MueLu_MapTransferFactory.hpp>
#include <MueLu_AggregationExportFactory.hpp>

#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListUtils.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include <MueLu_UseShortNames.hpp>

#include <MueLu_EpetraOperator.hpp> // Aztec interface

#include "muelu/muelu_ContactTransferFactory_decl.hpp"
#include "muelu/muelu_ContactASlaveDofFilterFactory_decl.hpp"
#include "muelu/MueLu_MyTrilinosSmoother_decl.hpp"

#include "solver_muelucontactpreconditioner.H"

// Read a parameter value from a parameter list and store it into a variable named 'varName'
#define MUELU_READ_PARAM(paramList, paramStr, varType, defaultValue, varName) \
  varType varName = defaultValue; if (paramList.isParameter(paramStr)) varName = paramList.get<varType>(paramStr);

// Read a parameter value from a paraeter list and copy it into a new parameter list (with another parameter name)
#define MUELU_COPY_PARAM(paramList, paramStr, varType, defaultValue, outParamList, outParamStr) \
  if (paramList.isParameter(paramStr))                                  \
    outParamList.set<varType>(outParamStr, paramList.get<varType>(paramStr)); \
  else outParamList.set<varType>(outParamStr, defaultValue);

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuContactPreconditioner::MueLuContactPreconditioner( FILE * outfile, Teuchos::ParameterList & mllist )
  : PreconditionerType( outfile ),
    mllist_( mllist )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuContactPreconditioner::Setup( bool create,
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
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO,LMO> >   mueluOp = Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO>(mueluA));

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
Teuchos::RCP<Hierarchy> LINALG::SOLVER::MueLuContactPreconditioner::SetupHierarchy(
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
  Teuchos::RCP<Epetra_Map> epMasterDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epSlaveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epActiveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epInnerDofMap = Teuchos::null;
  if(params.isSublist("Linear System properties")) {
    const Teuchos::ParameterList & linSystemProps = params.sublist("Linear System properties");
    //linSystemProps.set<Teuchos::RCP<Map> >("non diagonal-dominant row map",nonDiagMap);
    epMasterDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map> > ("contact masterDofMap");
    epSlaveDofMap =  linSystemProps.get<Teuchos::RCP<Epetra_Map> > ("contact slaveDofMap");
    epActiveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map> > ("contact activeDofMap");
    epInnerDofMap  = linSystemProps.get<Teuchos::RCP<Epetra_Map> > ("contact innerDofMap");
  }

  // build map extractor from different maps
  // note that the ordering (Master, Slave, Inner) is important to be the same overall the whole algorithm
  Teuchos::RCP<const Map> xfullmap = A->getRowMap(); // full map (MasterDofMap + SalveDofMap + InnerDofMap)
  //Teuchos::RCP<Xpetra::EpetraMap> xMasterDofMap  = Teuchos::rcp(new Xpetra::EpetraMap( epMasterDofMap ));
  Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epSlaveDofMap  ));
  //Teuchos::RCP<Xpetra::EpetraMap> xActiveDofMap  = Teuchos::rcp(new Xpetra::EpetraMap( epActiveDofMap ));
  //Teuchos::RCP<Xpetra::EpetraMap> xInnerDofMap   = Teuchos::rcp(new Xpetra::EpetraMap( epInnerDofMap  )); // TODO check me

  ///////////////////////////////////////////////////////////

  // calculate number of DOFs in xSlaveDofMap to extend maxCoarseSize
  Xpetra::global_size_t nSlaveDofs = xSlaveDofMap->getGlobalNumElements();

  ///////////////////////////////////////////////////////////

  // fill hierarchy
  Teuchos::RCP<Hierarchy> hierarchy = Teuchos::rcp(new Hierarchy(A));
  hierarchy->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(eVerbLevel));
  hierarchy->SetMaxCoarseSize(Teuchos::as<Xpetra::global_size_t>(maxCoarseSize)+nSlaveDofs);
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

 //Finest->Set("ActiveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO,GO,Node> >(xActiveDofMap));  // set map with active dofs
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
  Teuchos::RCP<UncoupledAggregationFactory> UCAggFact = Teuchos::rcp(new UncoupledAggregationFactory(dropFact));
  UCAggFact->SetMinNodesPerAggregate(minPerAgg);
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  UCAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  //UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);

  UCAggFact->SetOnePtMapName("SlaveDofMap", MueLu::NoFactory::getRCP());
  //UCAggFact->SetSmallAggMapName("SlaveDofMap", MueLu::NoFactory::getRCP());

  Teuchos::RCP<PFactory> PFact;
  Teuchos::RCP<TwoLevelFactoryBase> RFact;

  Teuchos::RCP<PFactory> PtentFact = Teuchos::rcp(new TentativePFactory());

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
    RFact = Teuchos::rcp( new TransPFactory() );
  } else {
    // Petrov Galerkin PG-AMG smoothed aggregation (energy minimization in ML)
    PFact  = Teuchos::rcp( new PgPFactory(/*PtentFact,slaveDcAFact*/) ); // use slaveDcAFact for prolongator smoothing
    PFact->SetFactory("P",PtentFact);
    PFact->SetFactory("A",slaveDcAFact);
    RFact  = Teuchos::rcp( new GenericRFactory() );
  }

  // define nullspace factory AFTER tentative PFactory (that generates the nullspace for the coarser levels)
  // use same nullspace factory for all multigrid levels
  // therefor we have to create one instance of NullspaceFactory and use it
  // for all FactoryManager objects (note: here, we have one FactoryManager object per level)
  Teuchos::RCP<NullspaceFactory> nspFact = Teuchos::rcp(new NullspaceFactory("Nullspace"/*,PtentFact*/));
  nspFact->SetFactory("Nullspace", PtentFact);

  // RAP factory with inter-level transfer of segregation block information (map extractor)
  Teuchos::RCP<RAPFactory> AcFact = Teuchos::rcp( new RAPFactory(/*PFact, RFact*/) );
  AcFact->SetFactory("P",PFact);
  AcFact->SetFactory("R",RFact);
  //AcFact->setVerbLevel(Teuchos::VERB_HIGH);
  AcFact->SetRepairZeroDiagonal(true); // repair zero diagonal entries in Ac, that are resulting from Ptent with nullspacedim > ndofspernode

  // write out aggregates
  Teuchos::RCP<MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > aggExpFact = Teuchos::rcp(new MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>());
  aggExpFact->SetParameter("Output filename",Teuchos::ParameterEntry(std::string("aggs_level%LEVELID_proc%PROCID.out")));
  aggExpFact->SetFactory("Aggregates",UCAggFact);
  aggExpFact->SetFactory("DofsPerNode",dropFact);
  //Input(fineLevel, "Aggregates");         //< factory which created aggregates
  //Input(fineLevel, "DofsPerNode");        //< CoalesceAndDropFactory (needed for DofsPerNode variable)
  //Input(fineLevel, "UnAmalgamationInfo"); //< AmalgamationFactory (needed for UnAmalgamationInfo variable)

  AcFact->AddTransferFactory(aggExpFact);

  // transfer maps to coarser grids
  //Teuchos::RCP<MueLu::MapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > cmTransFact3 = Teuchos::rcp(new MueLu::MapTransferFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", PtentFact, MueLu::NoFactory::getRCP()));
  Teuchos::RCP<MapTransferFactory> cmTransFact3 = Teuchos::rcp(new MapTransferFactory("SlaveDofMap", MueLu::NoFactory::getRCP()));
  cmTransFact3->SetFactory("P", PtentFact);
  AcFact->AddTransferFactory(cmTransFact3);

  ///////////////////////////////////////////////////////////////////////
  // setup coarse level smoothers/solvers
  ///////////////////////////////////////////////////////////////////////

  // coarse level smoother/solver
  //Teuchos::RCP<SmootherFactory> coarsestSmooFact;
  //coarsestSmooFact = MueLu::MLParameterListInterpreter<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetSmootherFactory(params, "coarse"); 
  Teuchos::RCP<SmootherFactory> coarsestSmooFact;
  coarsestSmooFact = GetContactCoarsestSolverFactory(params);

  ///////////////////////////////////////////////////////////////////////
  // prepare factory managers
  ///////////////////////////////////////////////////////////////////////

  bool bIsLastLevel = false;
  std::vector<Teuchos::RCP<FactoryManager> > vecManager(maxLevels);
  for(int i=0; i < maxLevels; i++) {
    //params.set("smoother: pre or post","pre"); // only pre-smoothing
    Teuchos::ParameterList pp(params);
    //pp.set("smoother: pre or post","pre");
    
    vecManager[i] = Teuchos::rcp(new FactoryManager());

    // fine/intermedium level smoother
    Teuchos::RCP<SmootherFactory> SmooFactFine = GetContactSmootherFactory(pp, i, slaveDcAFact);

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
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > LINALG::SOLVER::MueLuContactPreconditioner::GetContactSmootherFactory(const Teuchos::ParameterList & paramList, int level, const Teuchos::RCP<FactoryBase> & AFact) {

  char levelchar[11];
  sprintf(levelchar,"(level %d)",level);
  std::string levelstr(levelchar);

  if(paramList.isSublist("smoother: list " + levelstr)==false)
    return Teuchos::null;
  TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist("smoother: list " + levelstr)==false, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: no ML smoother parameter list for level. error.");

  std::string type = paramList.sublist("smoother: list " + levelstr).get<std::string>("smoother: type");
  TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: no ML smoother type for level. error.");

  const Teuchos::ParameterList smolevelsublist = paramList.sublist("smoother: list " + levelstr);
  //std::cout << "smoother: list " << levelstr << std::endl;
  //std::cout << smolevelsublist << std::endl;

  Teuchos::RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  Teuchos::RCP<SmootherFactory> SmooFact;

  if(type == "Jacobi") {
    if(smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", smolevelsublist.get<int>("smoother: sweeps"));
    if(smolevelsublist.get<double>("smoother: damping factor"))
      ifpackList.set("relaxation: damping factor", smolevelsublist.get<double>("smoother: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Jacobi");
    smooProto = Teuchos::rcp( new MueLu::MyTrilinosSmoother<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact) );
  } else if(type == "Gauss-Seidel") {
    if(smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", smolevelsublist.get<int>("smoother: sweeps"));
    if(smolevelsublist.get<double>("smoother: damping factor"))
      ifpackList.set("relaxation: damping factor", smolevelsublist.get<double>("smoother: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Gauss-Seidel");
    smooProto = Teuchos::rcp( new MueLu::MyTrilinosSmoother<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact) );
  } else if (type == "symmetric Gauss-Seidel") {
    if(smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", smolevelsublist.get<int>("smoother: sweeps"));
    if(smolevelsublist.get<double>("smoother: damping factor"))
      ifpackList.set("relaxation: damping factor", smolevelsublist.get<double>("smoother: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto = Teuchos::rcp( new MueLu::MyTrilinosSmoother<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact) );
    //std::cout << "built symm GS: " << smooProto << std::endl;
  } else if (type == "Chebyshev") {
    ifpackType = "CHEBYSHEV";
    if(smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set("chebyshev: degree", smolevelsublist.get<int>("smoother: sweeps"));
    smooProto = Teuchos::rcp( new MueLu::MyTrilinosSmoother<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact) );
    // TODO what about the other parameters
  } else if(type == "IFPACK") {
#ifdef HAVE_MUELU_IFPACK
    // TODO change to TrilinosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
    ifpackType = paramList.sublist("smoother: list " + levelstr).get<std::string>("smoother: ifpack type");
    if(ifpackType == "ILU") {
      ifpackList.set<int>("fact: level-of-fill", (int)smolevelsublist.get<double>("smoother: ifpack level-of-fill"));
      ifpackList.set("partitioner: overlap", smolevelsublist.get<int>("smoother: ifpack overlap"));
      int overlap = smolevelsublist.get<int>("smoother: ifpack overlap");
      //smooProto = MueLu::GetIfpackSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(ifpackType, ifpackList,smolevelsublist.get<int>("smoother: ifpack overlap"),AFact);
      smooProto = Teuchos::rcp( new MueLu::MyTrilinosSmoother<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>("SlaveDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, overlap, AFact) );
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: unknown ML smoother type " + type + " (IFPACK) not supported by MueLu. Only ILU is supported.");
#else // HAVE_MUELU_IFPACK
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: MueLu compiled without Ifpack support");
#endif // HAVE_MUELU_IFPACK
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: unknown ML smoother type " + type + " not supported by MueLu.");
  }

  // create smoother factory
  SmooFact = Teuchos::rcp( new SmootherFactory(smooProto) );

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if(smolevelsublist.isParameter("smoother: pre or post")) preorpost = smolevelsublist.get<std::string>("smoother: pre or post");

  if (preorpost == "pre") {
    SmooFact->SetSmootherPrototypes(smooProto, Teuchos::null);
  } else if(preorpost == "post") {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smooProto);
  }

  return SmooFact;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > LINALG::SOLVER::MueLuContactPreconditioner::GetContactCoarsestSolverFactory(const Teuchos::ParameterList & paramList, const Teuchos::RCP<FactoryBase> & AFact) {

  std::string type = ""; // use default defined by AmesosSmoother or Amesos2Smoother

  if(paramList.isParameter("coarse: type")) type = paramList.get<std::string>("coarse: type");

  Teuchos::RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  Teuchos::RCP<SmootherFactory> SmooFact;

  if(type == "Jacobi") {
    if(paramList.isParameter("coarse: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", paramList.get<int>("coarse: sweeps"));
    else ifpackList.set<int>("relaxation: sweeps", 1);
    if(paramList.isParameter("coarse: damping factor"))
      ifpackList.set("relaxation: damping factor", paramList.get<double>("coarse: damping factor"));
    else ifpackList.set("relaxation: damping factor", 1.0);
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Jacobi");
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList, 0, AFact) );
  } else if(type == "Gauss-Seidel") {
    if(paramList.isParameter("coarse: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", paramList.get<int>("coarse: sweeps"));
    else ifpackList.set<int>("relaxation: sweeps", 1);
    if(paramList.isParameter("coarse: damping factor"))
      ifpackList.set("relaxation: damping factor", paramList.get<double>("coarse: damping factor"));
    else ifpackList.set("relaxation: damping factor", 1.0);
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Gauss-Seidel");
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList, 0, AFact) );
  } else if (type == "symmetric Gauss-Seidel") {
    if(paramList.isParameter("coarse: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", paramList.get<int>("coarse: sweeps"));
    else ifpackList.set<int>("relaxation: sweeps", 1);
    if(paramList.isParameter("coarse: damping factor"))
      ifpackList.set("relaxation: damping factor", paramList.get<double>("coarse: damping factor"));
    else ifpackList.set("relaxation: damping factor", 1.0);
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList, 0, AFact) );
  } else if (type == "Chebyshev") {
    ifpackType = "CHEBYSHEV";
    if(paramList.isParameter("coarse: sweeps"))
      ifpackList.set("chebyshev: degree", paramList.get<int>("coarse: sweeps"));
    if(paramList.isParameter("coarse: Chebyshev alpha"))
      ifpackList.set("chebyshev: alpha", paramList.get<double>("coarse: Chebyshev alpha"));
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList, 0, AFact) );
  } else if(type == "IFPACK") {
#ifdef HAVE_MUELU_IFPACK
    // TODO change to TrilinosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
   /* ifpackType = paramList.get<std::string>("coarse: ifpack type");
    if(ifpackType == "ILU") {
      ifpackList.set<int>("fact: level-of-fill", (int)paramList.get<double>("coarse: ifpack level-of-fill"));
      ifpackList.set("partitioner: overlap", paramList.get<int>("coarse: ifpack overlap"));
      smooProto = MueLu::GetIfpackSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(ifpackType, ifpackList, paramList.get<int>("coarse: ifpack overlap"), AFact);
    }
    else*/
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: unknown ML smoother type " + type + " (IFPACK) not supported by MueLu. Only ILU is supported.");
#else // HAVE_MUELU_IFPACK
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: MueLu compiled without Ifpack support");
#endif // HAVE_MUELU_IFPACK
  } else if(type == "Amesos-Superlu") {
    smooProto = Teuchos::rcp( new DirectSolver("Superlu",Teuchos::ParameterList()) );
  } else if(type == "Amesos-Superludist") {
    smooProto = Teuchos::rcp( new DirectSolver("Superludist",Teuchos::ParameterList()) );
  } else if(type == "Amesos-KLU") {
    smooProto = Teuchos::rcp( new DirectSolver("Klu",Teuchos::ParameterList()) );
  } else if(type == "Amesos-UMFPACK") {
    smooProto = Teuchos::rcp( new DirectSolver("Umfpack",Teuchos::ParameterList()) );
  } else if(type == "") {
    smooProto = Teuchos::rcp( new DirectSolver("",Teuchos::ParameterList()) );
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: unknown coarsest solver type. '" << type << "' not supported by MueLu.");
  }

  // set AFactory
  smooProto->SetFactory("A", AFact);

  // create smoother factory
  TEUCHOS_TEST_FOR_EXCEPTION(smooProto == Teuchos::null, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter: no smoother prototype. fatal error.");
  SmooFact = rcp( new SmootherFactory(smooProto) );

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if(paramList.isParameter("coarse: pre or post")) preorpost = paramList.get<std::string>("coarse: pre or post");

  if (preorpost == "pre") {
    SmooFact->SetSmootherPrototypes(smooProto, Teuchos::null);
  } else if(preorpost == "post") {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smooProto);
  }

  return SmooFact;
}

#endif // HAVE_EXPERIMENTAL_MueLu
#endif

