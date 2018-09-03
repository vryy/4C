/*!----------------------------------------------------------------------
\file solver_muelucontactpreconditioner_old.cpp

\brief Implementation

\level 1

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/
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

#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_UCAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_PgPFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_NullspaceFactory.hpp>
//#include <MueLu_SegregationAFilterFactory.hpp>
#include <MueLu_SegregationATransferFactory.hpp>  // TODO remove me
#include <MueLu_Aggregates.hpp>

#include <MueLu_AggregationExportFactory.hpp>

#include <MueLu_MLParameterListInterpreter_decl.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include <MueLu_EpetraOperator.hpp>  // Aztec interface

#include "muelu_ContactAFilterFactory_decl.hpp"
#include "muelu_ContactTransferFactory_decl.hpp"

#include "solver_muelucontactpreconditioner.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuContactPreconditioner::MueLuContactPreconditioner(
    FILE* outfile, Teuchos::ParameterList& mllist)
    : PreconditionerType(outfile), mllist_(mllist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuContactPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  if (create)
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
    if (A == NULL) dserror("CrsMatrix expected");

    // free old matrix first
    P_ = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

    // see whether we use standard ml or our own mlapi operator
    // const bool domuelupreconditioner = mllist_.get<bool>("LINALG::MueLu_Preconditioner",false);

    // wrap Epetra_CrsMatrix to Xpetra::Operator for use in MueLu
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO, LMO>> mueluA =
        Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
    Teuchos::RCP<Xpetra::Operator<SC, LO, GO, NO, LMO>> mueluOp =
        Teuchos::rcp(new Xpetra::CrsOperator<SC, LO, GO, NO, LMO>(mueluA));

    // prepare nullspace vector for MueLu
    int numdf = mllist_.get<int>("PDE equations", -1);
    int dimns = mllist_.get<int>("null space: dimension", -1);
    if (dimns == -1 || numdf == -1) dserror("Error: PDE equations or null space dimension wrong.");
    Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> rowMap = mueluA->getRowMap();

    Teuchos::RCP<MultiVector> nspVector =
        Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
            rowMap, dimns, true);
    Teuchos::RCP<std::vector<double>> nsdata =
        mllist_.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);

    for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
    {
      Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
      const size_t myLength = nspVector->getLocalLength();
      for (size_t j = 0; j < myLength; j++)
      {
        nspVectori[j] = (*nsdata)[i * myLength + j];
      }
    }

    // remove unsupported flags
    mllist_.remove("aggregation: threshold", false);  // no support for aggregation: threshold TODO

    // Setup MueLu Hierarchy
    // Teuchos::RCP<Hierarchy> H = MLInterpreter::Setup(mllist_, mueluOp, nspVector);
    Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H =
        SetupHierarchy(mllist_, mueluOp, nspVector);

    // set preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Hierarchy> LINALG::SOLVER::MueLuContactPreconditioner::SetupHierarchy(
    const Teuchos::ParameterList& params, const Teuchos::RCP<Operator>& A,
    const Teuchos::RCP<MultiVector> nsp)
{
  //#include "MueLu_UseShortNames.hpp" // TODO don't know why this is needed here

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // read in common parameters
  int maxLevels = 10;       // multigrid prameters
  int verbosityLevel = 10;  // verbosity level
  int maxCoarseSize = 50;
  int nDofsPerNode = 1;        // coalesce and drop parameters
  double agg_threshold = 0.0;  // aggregation parameters
  double agg_damping = 4 / 3;
  int agg_smoothingsweeps = 1;
  int minPerAgg = 3;  // optimal for 2d
  int maxNbrAlreadySelected = 0;
  std::string agg_type = "Uncoupled";
  bool bEnergyMinimization = false;  // PGAMG
  if (params.isParameter("max levels")) maxLevels = params.get<int>("max levels");
  if (params.isParameter("ML output")) verbosityLevel = params.get<int>("ML output");
  if (params.isParameter("coarse: max size")) maxCoarseSize = params.get<int>("coarse: max size");
  if (params.isParameter("PDE equations")) nDofsPerNode = params.get<int>("PDE equations");
  if (params.isParameter("aggregation: threshold"))
    agg_threshold = params.get<double>("aggregation: threshold");
  if (params.isParameter("aggregation: damping factor"))
    agg_damping = params.get<double>("aggregation: damping factor");
  if (params.isParameter("aggregation: smoothing sweeps"))
    agg_smoothingsweeps = params.get<int>("aggregation: smoothing sweeps");
  if (params.isParameter("aggregation: type"))
    agg_type = params.get<std::string>("aggregation: type");
  if (params.isParameter("aggregation: nodes per aggregate"))
    minPerAgg = params.get<int>("aggregation: nodes per aggregate");
  if (params.isParameter("energy minimization: enable"))
    bEnergyMinimization = params.get<bool>("energy minimization: enable");

  // set DofsPerNode in A operator
  A->SetFixedBlockSize(nDofsPerNode);

  // translate verbosity parameter
  Teuchos::EVerbosityLevel eVerbLevel = Teuchos::VERB_NONE;
  if (verbosityLevel == 0) eVerbLevel = Teuchos::VERB_NONE;
  if (verbosityLevel > 0) eVerbLevel = Teuchos::VERB_LOW;
  if (verbosityLevel > 4) eVerbLevel = Teuchos::VERB_MEDIUM;
  if (verbosityLevel > 7) eVerbLevel = Teuchos::VERB_HIGH;
  if (verbosityLevel > 9) eVerbLevel = Teuchos::VERB_EXTREME;

  // extract additional maps from parameter list
  // these maps are provided by the STR::TimInt::PrepareContactMeshtying routine, that
  // has access to the contact manager class
  Teuchos::RCP<Epetra_Map> epMasterDofMap = params.get<Teuchos::RCP<Epetra_Map>>(
      "LINALG::SOLVER::MueLu_ContactPreconditioner::MasterDofMap");
  Teuchos::RCP<Epetra_Map> epSlaveDofMap = params.get<Teuchos::RCP<Epetra_Map>>(
      "LINALG::SOLVER::MueLu_ContactPreconditioner::SlaveDofMap");
  // Teuchos::RCP<Epetra_Map> epInnerDofMap  = params.get<Teuchos::RCP<Epetra_Map>
  // >("LINALG::SOLVER::MueLu_ContactPreconditioner::InnerDofMap"); // TODO check me

  // build map extractor from different maps
  // note that the ordering (Master, Slave, Inner) is important to be the same overall the whole
  // algorithm
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> xfullmap =
      A->getRowMap();  // full map (MasterDofMap + SalveDofMap + InnerDofMap)
  Teuchos::RCP<Xpetra::EpetraMap> xMasterDofMap =
      Teuchos::rcp(new Xpetra::EpetraMap(epMasterDofMap));
  Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap = Teuchos::rcp(new Xpetra::EpetraMap(epSlaveDofMap));
  // Teuchos::RCP<Xpetra::EpetraMap> xInnerDofMap   = Teuchos::rcp(new Xpetra::EpetraMap(
  // epInnerDofMap  )); // TODO check me

  std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, Node>>> xmaps;
  xmaps.push_back(xMasterDofMap);
  xmaps.push_back(xSlaveDofMap);
  // xmaps.push_back(xInnerDofMap ); // TODO check me

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node>> map_extractor =
      Xpetra::MapExtractorFactory<Scalar, LO, GO>::Build(xfullmap, xmaps);

  // create factories

  // prepare (filtered) A Factory
  Teuchos::RCP<MueLu::SingleLevelFactoryBase> segAFact = Teuchos::rcp(
      new MueLu::ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(
          "A", NULL, map_extractor));

  // nullspace factory
  Teuchos::RCP<NullspaceFactory> nspFact = Teuchos::rcp(new NullspaceFactory());

  // Coalesce and drop factory with constant number of Dofs per freedom
  Teuchos::RCP<CoalesceDropFactory> dropFact =
      Teuchos::rcp(new CoalesceDropFactory(segAFact, nspFact));
  // dropFact->SetFixedBlockSize(nDofsPerNode);

  // aggregation factory
  Teuchos::RCP<MueLu::UCAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> UCAggFact =
      Teuchos::rcp(
          new MueLu::UCAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>(dropFact));
  // note: this class does not derive from VerboseObject. Therefore we cannot use GetOStream
  if (verbosityLevel > 3)
  {
    *out << "========================= Aggregate option summary  ========================="
         << std::endl;
    *out << "min Nodes per aggregate :               " << minPerAgg << std::endl;
    *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
    *out << "aggregate ordering :                    NATURAL" << std::endl;
    *out << "============================================================================="
         << std::endl;
  }
  UCAggFact->SetMinNodesPerAggregate(
      minPerAgg);  // TODO should increase if run anything other than 1D
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  UCAggFact->SetPhase3AggCreation(0.5);

  // transfer operators (PG-AMG)
  /*Teuchos::RCP<PFactory> PtentFact = Teuchos::rcp(new
  TentativePFactory(UCAggFact,nspFact,segAFact)); Teuchos::RCP<PFactory> PFact  = Teuchos::rcp( new
  PgPFactory(PtentFact,segAFact) ); Teuchos::RCP<RFactory> RFact  = Teuchos::rcp( new
  GenericRFactory(PFact) );*/
  Teuchos::RCP<MueLu::PFactory> PFact =
      Teuchos::rcp(new TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
          UCAggFact, nspFact, segAFact));
  Teuchos::RCP<MueLu::RFactory> RFact =
      Teuchos::rcp(new TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>(PFact));

  // RAP factory with inter-level transfer of segregation block information (map extractor)
  Teuchos::RCP<MueLu::RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AcFact =
      Teuchos::rcp(new MueLu::RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>(PFact, RFact));

  // write out aggregates
  Teuchos::RCP<
      MueLu::AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>>
      aggExpFact = Teuchos::rcp(new MueLu::AggregationExportFactory<Scalar, LocalOrdinal,
          GlobalOrdinal, Node, LocalMatOps>(
          "aggs_level%LEVELID_proc%PROCID.out", UCAggFact.get(), dropFact.get(), segAFact.get()));
  AcFact->AddTransferFactory(aggExpFact);

  Teuchos::RCP<
      MueLu::ContactTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>>
      cTransFact = Teuchos::rcp(
          new MueLu::ContactTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(
              PFact));
  AcFact->AddTransferFactory(cTransFact);

  // setup smoothers
  Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> coarsestSmooFact;
  coarsestSmooFact = MLParameterListInterpreter::GetCoarsestSolverFactory(params);

  ///////////////////////////////////////////////////

  // fill hierarchy
  Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> hierarchy =
      Teuchos::rcp(new MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A));
  hierarchy->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(eVerbLevel));
  hierarchy->SetMaxCoarseSize(Teuchos::as<Xpetra::global_size_t>(maxCoarseSize));

  ///////////////////////////////////////////////////////////

  // set fine level nullspace
  // use given fine level null space or extract pre-computed nullspace from ML parameter list
  Teuchos::RCP<MueLu::Level> Finest = hierarchy->GetLevel();  // get finest level
  if (nsp != Teuchos::null)
  {
    Finest->Set("Nullspace", nsp);  // set user given null space
  }
  else
  {
    std::string type = "";
    if (params.isParameter("null space: type")) type = params.get<std::string>("null space: type");
    if (type != "pre-computed")
      dserror("MueLu::Interpreter: no valid nullspace (no pre-computed null space). error.");
    int dimns = -1;
    if (params.isParameter("null space: dimension"))
      dimns = params.get<int>("null space: dimension");
    if (dimns == -1) dserror("MueLu::Interpreter: no valid nullspace (nullspace dim = -1). error.");

    const Teuchos::RCP<const Map> rowMap = A->getRowMap();
    Teuchos::RCP<MultiVector> nspVector = MultiVectorFactory::Build(rowMap, dimns, true);
    double* nsdata = NULL;
    if (params.isParameter("null space: vectors"))
      nsdata = params.get<double*>("null space: vectors");
    if (nsdata == NULL) dserror("MueLu::Interpreter: no valid nullspace (nsdata = NULL). error.");

    for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
    {
      Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
      const size_t myLength = nspVector->getLocalLength();
      for (size_t j = 0; j < myLength; j++)
      {
        nspVectori[j] = nsdata[i * myLength + j];
      }
    }

    Finest->Set("Nullspace", nspVector);  // set user given null space
  }

  Finest->Set("A", A);
  Finest->Keep("Aggregates", UCAggFact.get());
  Finest->Keep("A", segAFact.get());  // TODO: keep segAfact A for export of aggregates!
  ////////////////////////////////////

  // prepare factory managers

  bool bIsLastLevel = false;
  std::vector<Teuchos::RCP<MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>
      vecManager(maxLevels);
  for (int i = 0; i < maxLevels; i++)
  {
    Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> SmooFactFine =
        MLParameterListInterpreter::GetSmootherFactory(params, i);

    vecManager[i] =
        Teuchos::rcp(new MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
    if (SmooFactFine != Teuchos::null)
      vecManager[i]->SetFactory("Smoother",
          SmooFactFine);  // Hierarchy.Setup uses TOPSmootherFactory, that only needs "Smoother"
    vecManager[i]->SetFactory("CoarseSolver", coarsestSmooFact);
    vecManager[i]->SetFactory("A", AcFact);  // same RAP factory
    vecManager[i]->SetFactory("P", PFact);   // same prolongator and restrictor factories
    vecManager[i]->SetFactory("R", RFact);   // same prolongator and restrictor factories
    vecManager[i]->SetFactory(
        "Nullspace", nspFact);  // use same nullspace factory throughout all multigrid levels
  }

  // use new Hierarchy::Setup routine
  if (maxLevels == 1)
  {
    bIsLastLevel = hierarchy->Setup(0, Teuchos::null, vecManager[0].ptr(), Teuchos::null);
  }
  else
  {
    bIsLastLevel = hierarchy->Setup(0, Teuchos::null, vecManager[0].ptr(),
        vecManager[1].ptr());  // true, false because first level
    for (int i = 1; i < maxLevels - 1; i++)
    {
      if (bIsLastLevel == true) break;
      bIsLastLevel = hierarchy->Setup(
          i, vecManager[i - 1].ptr(), vecManager[i].ptr(), vecManager[i + 1].ptr());
    }
    if (bIsLastLevel == false)
    {
      bIsLastLevel = hierarchy->Setup(maxLevels - 1, vecManager[maxLevels - 2].ptr(),
          vecManager[maxLevels - 1].ptr(), Teuchos::null);
    }
  }


  return hierarchy;
}



#endif
