/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de

Created on: Sep 7, 2012
Author: wiesner
*----------------------------------------------------------------------*/
#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)

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
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_NullspaceFactory.hpp>
#include <MueLu_Aggregates.hpp>
#include <MueLu_MapTransferFactory.hpp>
#include <MueLu_AggregationExportFactory.hpp>
#include <MueLu_IfpackSmoother.hpp>

#include <MueLu_MLParameterListInterpreter.hpp>

#ifdef HAVE_MUELU_ISORROPIA

#include "MueLu_IsorropiaInterface.hpp"
#include "MueLu_RepartitionInterface.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceMapFactory.hpp"
#endif



// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include <MueLu_EpetraOperator.hpp>  // Aztec interface

#include "muelu/MueLu_ContactAFilterFactory_decl.hpp"
#include "muelu/MueLu_ContactTransferFactory_decl.hpp"
#include "muelu/MueLu_MyTrilinosSmoother_decl.hpp"
#include "muelu/MueLu_IterationAFactory_decl.hpp"
#include "muelu/MueLu_SelectiveSaPFactory_decl.hpp"

#include "solver_muelucontactpreconditioner2.H"

typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;

typedef MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> Hierarchy;
typedef MueLu::Factory Factory;
typedef MueLu::CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CoalesceDropFactory;
typedef MueLu::UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>
    UncoupledAggregationFactory;
typedef MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> SmootherFactory;
typedef MueLu::TentativePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> TentativePFactory;
typedef MueLu::TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> TransPFactory;
typedef MueLu::NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> NullspaceFactory;
typedef MueLu::GenericRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> GenericRFactory;
typedef MueLu::MapTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapTransferFactory;
typedef MueLu::RebalanceAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> RebalanceAcFactory;
typedef MueLu::AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> AmalgamationFactory;
typedef MueLu::PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> PgPFactory;
typedef MueLu::RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> RAPFactory;
typedef MueLu::RepartitionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> RepartitionFactory;
typedef MueLu::RebalanceTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    RebalanceTransferFactory;
typedef MueLu::RebalanceMapFactory<LocalOrdinal, GlobalOrdinal, Node> RebalanceMapFactory;

typedef Scalar SC;
typedef LocalOrdinal LO;
typedef GlobalOrdinal GO;
typedef Node NO;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuContactPreconditioner2::MueLuContactPreconditioner2(
    FILE* outfile, Teuchos::ParameterList& mllist)
    : PreconditionerType(outfile), mllist_(mllist)
{
  singleNodeAFact_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuContactPreconditioner2::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  // handle reuse strategy
  std::string reuseStrategy = "nothing";
  if (mllist_.isParameter("muelu reuse: strategy"))
    reuseStrategy = mllist_.get<std::string>("muelu reuse: strategy");

  // reuseStrategy == "nothing" -> force rebuilding the preconditioner
  if (reuseStrategy == "nothing" && create == false) create = true;

  if (create)
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
    if (A == NULL) dserror("CrsMatrix expected");

    // free old matrix first
    P_ = Teuchos::null;
    H_ = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

    // wrap Epetra_CrsMatrix to Xpetra::Operator for use in MueLu
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> mueluA =
        Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
    Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> mueluOp =
        Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(mueluA));

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

    H_ = Teuchos::rcp(new Hierarchy());
    H_ = InitializeHierarchy(H_, mllist_, mueluOp, nspVector);
    H_ = SetupFactories(H_, mllist_);
    H_ = SetupSmoothers(H_, mllist_);

    // Setup MueLu Hierarchy
    H_ = SetupHierarchy(H_, mllist_, mueluOp, nspVector);

    // set preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H_));
  }
  // only rebuild parts of preconditioner
  // reuse nonsmoothed transfer operators and Importer objects in case of active rebalancing
  else if (create == false && reuseStrategy == "Ptent")
  {
    H_->ExpertClear();

    /*Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); Teuchos::RCP<Level> coarseLevel0 =
    H_->GetLevel(0); coarseLevel0->print(*out);*/

    // wrap Epetra_CrsMatrix to Xpetra::Operator for use in MueLu
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> mueluA =
        Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
    Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> mueluOp =
        Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(mueluA));

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

    H_ = InitializeHierarchy(H_, mllist_, mueluOp, nspVector);

    H_ = SetupSmoothers(H_, mllist_);

    // Setup MueLu Hierarchy
    H_ = SetupHierarchy(H_, mllist_, mueluOp, nspVector);

    // set preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H_));
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>
LINALG::SOLVER::MueLuContactPreconditioner2::InitializeHierarchy(
    const Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>& h, const Teuchos::ParameterList& params,
    const Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>>& A,
    const Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nsp)
{
  int verbosityLevel = 10;  // verbosity level
  int nDofsPerNode = 1;     // coalesce and drop parameters
  int maxCoarseSize = 50;
  // bool bSegregateAggregates = true; //false; // TODO fix me: only for tests with repartitioning!
  // true; // segregate aggregates (to enforce non-overlapping aggregates between master and slave
  // side)

  if (params.isParameter("ML output")) verbosityLevel = params.get<int>("ML output");
  if (params.isParameter("PDE equations")) nDofsPerNode = params.get<int>("PDE equations");
  if (params.isParameter("coarse: max size")) maxCoarseSize = params.get<int>("coarse: max size");

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
  Teuchos::RCP<Epetra_Map> epMasterDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epSlaveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epActiveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epInnerDofMap = Teuchos::null;
  Teuchos::RCP<Map> xSingleNodeAggMap = Teuchos::null;
  Teuchos::RCP<Map> xNearZeroDiagMap = Teuchos::null;
  if (params.isSublist("Linear System properties"))
  {
    const Teuchos::ParameterList& linSystemProps = params.sublist("Linear System properties");
    // extract information provided by solver (e.g. PermutedAztecSolver)
    epMasterDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact masterDofMap");
    epSlaveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact slaveDofMap");
    epActiveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact activeDofMap");
    epInnerDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact innerDofMap");
    if (linSystemProps.isParameter("non diagonal-dominant row map"))
      xSingleNodeAggMap = linSystemProps.get<Teuchos::RCP<Map>>("non diagonal-dominant row map");
    if (linSystemProps.isParameter("near-zero diagonal row map"))
      xNearZeroDiagMap = linSystemProps.get<Teuchos::RCP<Map>>("near-zero diagonal row map");
    /*if(linSystemProps.isParameter("ProblemType") && linSystemProps.get<std::string>("ProblemType")
      == "meshtying") bSegregateAggregates = false;*/
  }

  // transform Epetra maps to Xpetra maps (if necessary)
  Teuchos::RCP<const Map> xfullmap =
      A->getRowMap();  // full map (MasterDofMap + SalveDofMap + InnerDofMap)
  Teuchos::RCP<Xpetra::EpetraMap> xMasterDofMap =
      Teuchos::rcp(new Xpetra::EpetraMap(epMasterDofMap));
  Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap = Teuchos::rcp(new Xpetra::EpetraMap(epSlaveDofMap));

  ///////////////////////////////////////////////////////////

  // fill hierarchy
  // Teuchos::RCP<Hierarchy> hierarchy = Teuchos::rcp(new Hierarchy(A));
  h->setlib(Xpetra::UseEpetra);
  h->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(eVerbLevel));
  h->SetMaxCoarseSize(Teuchos::as<Xpetra::global_size_t>(maxCoarseSize) /*+nSlaveDofs*/);

  ///////////////////////////////////////////////////////////

  // set fine level nullspace
  // use given fine level null space or extract pre-computed nullspace from ML parameter list
  Teuchos::RCP<MueLu::Level> Finest = h->GetLevel(0);  // get finest level
  Finest->setlib(Xpetra::UseEpetra);
  Finest->Set("A", A);  // not necessary


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

  ///////////////////////////////////////////////////////////////////////
  // declare "SingleNodeAggDofMap" on finest level
  // this map marks the problematic rows of the input matrix
  // these problematic rows are skipped within the multigrid smoothers
  // -> MyTrilinosSmoother, IterationAFactory
  if (xSingleNodeAggMap != Teuchos::null)
    Finest->Set("SingleNodeAggDofMap",
        Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, Node>>(xSingleNodeAggMap));

  ///////////////////////////////////////////////////////////////////////
  // declare "SlaveDofMap" on finest level
  // this map marks the slave DOFs (ignoring permutations) which shall
  // be excluded from transfer operator smoothing
  // -> SelectiveSaPFactory
  if (xSlaveDofMap != Teuchos::null)
    Finest->Set(
        "SlaveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, Node>>(xSlaveDofMap));

  ///////////////////////////////////////////////////////////////////////
  // declare "MasterDofMap" on finest level
  // this map marks the master DOFs (ignoring permutations)
  // needed together with SlaveDofMap to define a segregated matrix, used for
  // building aggregates
  // -> ContactAFilterFactory
  if (xMasterDofMap != Teuchos::null)
    Finest->Set(
        "MasterDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, Node>>(xMasterDofMap));

  ///////////////////////////////////////////////////////////////////////
  // declare "NearZeroDiagMap" on finest level
  // this map marks rows of A with a near zero diagonal entry (usually the
  // tolerance is 1e-5). This is a sign of a badly behaving matrix. Usually
  // one should not used smoothed aggregation for these type of problems but
  // rather try plain aggregation (PA-AMG)
  // -> SelectiveSaPFactory
  if (xNearZeroDiagMap != Teuchos::null)
    Finest->Set("NearZeroDiagMap",
        Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, Node>>(xNearZeroDiagMap));

  // now, the finest level contain the following variables
  // - A
  // - Nullspace
  // - NearZeroDiagMap (if available)
  // - MasterDofMap (if available)
  // - SlaveDofMap (if available)
  // - SingleNodeAggDofMap (if available)

  return h;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>
LINALG::SOLVER::MueLuContactPreconditioner2::SetupFactories(
    const Teuchos::RCP<Hierarchy>& h, const Teuchos::ParameterList& params)
{
  // read in common parameters
  int maxLevels = 10;  // multigrid prameters
  // double agg_threshold = 0.0;   // aggregation parameters
  double agg_damping = 4 / 3;
  int minPerAgg = 3;   // optimal for 2d
  int maxPerAgg = 27;  // optimal for 3d
  int maxNbrAlreadySelected = 0;
  std::string agg_type = "Uncoupled";
  std::string reuseStrategy = "none";
  bool bSegregateAggregates =
      true;  // false; // TODO fix me: only for tests with repartitioning! true; // segregate
             // aggregates (to enforce non-overlapping aggregates between master and slave side)

#ifdef HAVE_MUELU_ISORROPIA
  bool bDoRepartition = false;
  double optNnzImbalance = 1.3;
  int optMinRowsPerProc = 3000;
#endif

  if (params.isParameter("max levels")) maxLevels = params.get<int>("max levels");
  if (params.isParameter("aggregation: damping factor"))
    agg_damping = params.get<double>("aggregation: damping factor");
  // if(params.isParameter("aggregation: smoothing sweeps"))   agg_smoothingsweeps = params.get<int>
  // ("aggregation: smoothing sweeps");
  if (params.isParameter("aggregation: type"))
    agg_type = params.get<std::string>("aggregation: type");
  if (params.isParameter("aggregation: min nodes per aggregate"))
    minPerAgg = params.get<int>("aggregation: min nodes per aggregate");
  if (params.isParameter("aggregation: nodes per aggregate"))
    maxPerAgg = params.get<int>("aggregation: nodes per aggregate");
  if (params.isParameter("muelu reuse: strategy"))
    reuseStrategy = params.get<std::string>("muelu reuse: strategy");
#ifdef HAVE_MUELU_ISORROPIA
  if (params.isParameter("muelu repartition: enable"))
  {
    if (params.get<int>("muelu repartition: enable") == 1) bDoRepartition = true;
  }
  if (params.isParameter("muelu repartition: max min ratio"))
    optNnzImbalance = params.get<double>("muelu repartition: max min ratio");
  if (params.isParameter("muelu repartition: min per proc"))
    optMinRowsPerProc = params.get<int>("muelu repartition: min per proc");
#else
  dserror("Isorropia has not been compiled with Trilinos. Repartitioning is not working.");
#endif

  // std::cout << "Reuse strategy: " << reuseStrategy << std::endl;

  // extract additional maps from parameter list
  // these maps are provided by the STR::TimInt::PrepareContactMeshtying routine, that
  // has access to the contact manager class
  Teuchos::RCP<Epetra_Map> epMasterDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epSlaveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epActiveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epInnerDofMap = Teuchos::null;
  Teuchos::RCP<Map> xSingleNodeAggMap = Teuchos::null;
  Teuchos::RCP<Map> xNearZeroDiagMap = Teuchos::null;
  if (params.isSublist("Linear System properties"))
  {
    const Teuchos::ParameterList& linSystemProps = params.sublist("Linear System properties");
    // extract information provided by solver (e.g. PermutedAztecSolver)
    epMasterDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact masterDofMap");
    epSlaveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact slaveDofMap");
    epActiveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact activeDofMap");
    epInnerDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact innerDofMap");
    if (linSystemProps.isParameter("non diagonal-dominant row map"))
      xSingleNodeAggMap = linSystemProps.get<Teuchos::RCP<Map>>("non diagonal-dominant row map");
    if (linSystemProps.isParameter("near-zero diagonal row map"))
      xNearZeroDiagMap = linSystemProps.get<Teuchos::RCP<Map>>("near-zero diagonal row map");
    if (linSystemProps.isParameter("ProblemType") &&
        linSystemProps.get<std::string>("ProblemType") == "meshtying")
      bSegregateAggregates = false;
  }

  // transform Epetra maps to Xpetra maps (if necessary)
  // Teuchos::RCP<const Map> xfullmap = A->getRowMap(); // full map (MasterDofMap + SalveDofMap +
  // InnerDofMap)
  Teuchos::RCP<Xpetra::EpetraMap> xMasterDofMap =
      Teuchos::rcp(new Xpetra::EpetraMap(epMasterDofMap));
  Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap = Teuchos::rcp(new Xpetra::EpetraMap(epSlaveDofMap));

  ///////////////////////////////////////////////////////////

  Teuchos::RCP<MueLu::Level> Finest = h->GetLevel(0);  // get finest level


  // for segregating aggregates (slave and master)
  Teuchos::RCP<Factory> segregatedAFact = Teuchos::null;
  if (bSegregateAggregates)
  {
    segregatedAFact =
        Teuchos::rcp(new MueLu::ContactAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
    segregatedAFact->SetParameter("Input matrix name", Teuchos::ParameterEntry(std::string("A")));
    segregatedAFact->SetParameter(
        "Map block 1 name", Teuchos::ParameterEntry(std::string("SlaveDofMap")));
    segregatedAFact->SetParameter(
        "Map block 2 name", Teuchos::ParameterEntry(std::string("MasterDofMap")));
    segregatedAFact->SetParameter(
        "Map block 1 factory", Teuchos::ParameterEntry(std::string("NoFactory")));
    segregatedAFact->SetParameter(
        "Map block 2 factory", Teuchos::ParameterEntry(std::string("NoFactory")));
  }

  // Coalesce and drop factory with constant number of Dofs per freedom
  // note: coalescing based on original matrix A
  Teuchos::RCP<CoalesceDropFactory> dropFact = Teuchos::rcp(new CoalesceDropFactory());
  if (bSegregateAggregates)
    dropFact->SetFactory(
        "A", segregatedAFact);  // if segregated aggregates are wished, set A factory here

  // aggregation factory
  Teuchos::RCP<UncoupledAggregationFactory> UCAggFact =
      Teuchos::rcp(new UncoupledAggregationFactory());
  UCAggFact->SetFactory("Graph", dropFact);
  UCAggFact->SetFactory("DofsPerNode", dropFact);
  UCAggFact->SetParameter(
      "aggregation: max selected neighbors", Teuchos::ParameterEntry(maxNbrAlreadySelected));
  UCAggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(minPerAgg));
  UCAggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(maxPerAgg));
  UCAggFact->SetParameter("aggregation: ordering", Teuchos::ParameterEntry(std::string("graph")));

  if (xSingleNodeAggMap != Teuchos::null)
  {  // declare single node aggregates
    UCAggFact->SetParameter(
        "aggregation: allow user-specified singletons", Teuchos::ParameterEntry(true));
    UCAggFact->SetParameter(
        "OnePt aggregate map name", Teuchos::ParameterEntry(std::string("SingleNodeAggDofMap")));
    UCAggFact->SetParameter(
        "OnePt aggregate map factory", Teuchos::ParameterEntry(std::string("NoFactory")));
  }

  Teuchos::RCP<MueLu::PFactory> PFact;
  Teuchos::RCP<MueLu::TwoLevelFactoryBase> RFact;

  Teuchos::RCP<MueLu::PFactory> PtentFact = Teuchos::rcp(new TentativePFactory());

  if (reuseStrategy == "Ptent") h->Keep("P", PtentFact.get());  // keep tentative P factory!!!!

  // choose either nonsmoothed transfer operators or
  // PG-AMG smoothed aggregation transfer operators
  // note:
  //  - SA-AMG is not working properly (probably due to problematic Dinv scaling with zeros on
  //  diagonal) TODO handling of zeros on diagonal in SaPFactory
  //  - PG-AMG has some special handling for zeros on diagonal (avoid NaNs)
  //    avoid local damping factory omega==1 -> oversmoothing, leads to zero rows in P
  //    use matrix A with artificial Dirichlet bcs for prolongator smoothing
  // if agg_damping == 0.0 -> PA-AMG else PG-AMG
  if (agg_damping == 0.0)
  {
    // tentative prolongation operator (PA-AMG)
    PFact = PtentFact;
    RFact = Teuchos::rcp(new TransPFactory());
  }
  else if (agg_damping > 0.0)
  {
    // smoothed aggregation (SA-AMG)
    PFact =
        Teuchos::rcp(new MueLu::SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());

    // feed SelectiveSAPFactory with information
    PFact->SetFactory("P", PtentFact);
    // use user-given damping parameter
    PFact->SetParameter("Damping factor", Teuchos::ParameterEntry(agg_damping));
    PFact->SetParameter("Damping strategy", Teuchos::ParameterEntry(std::string("User")));
    // only prolongator smoothing for transfer operator basis functions which
    // correspond to non-slave rows in (permuted) matrix A
    // We use the tentative prolongator to detect the corresponding prolongator basis functions for
    // given row gids. Note: this ignores the permutations in A. In case, the matrix A has been
    // permuted it can happen
    //       that problematic columns in Ptent are not corresponding to columns that belong to the
    //       with nonzero entries in slave rows. // TODO think more about this -> aggregation
    PFact->SetParameter("NonSmoothRowMapName", Teuchos::ParameterEntry(std::string("SlaveDofMap")));
    PFact->SetFactory("NonSmoothRowMapFactory", MueLu::NoFactory::getRCP());

    // provide diagnostics of diagonal entries of current matrix A
    // if the solver object detects some significantly small entries on diagonal the contact
    // preconditioner can decide to skip transfer operator smoothing to increase robustness
    PFact->SetParameter(
        "NearZeroDiagMapName", Teuchos::ParameterEntry(std::string("NearZeroDiagMap")));
    PFact->SetFactory("NearZeroDiagMapFactory", MueLu::NoFactory::getRCP());

    if (bSegregateAggregates)
      PFact->SetFactory("A", segregatedAFact);  // make sure that prolongator smoothing does not
                                                // disturb segregation of transfer operators

    RFact = Teuchos::rcp(new GenericRFactory());
  }
  else
  {
    // Petrov Galerkin PG-AMG smoothed aggregation (energy minimization in ML)
    PFact = Teuchos::rcp(new PgPFactory());
    PFact->SetFactory("P", PtentFact);
    if (bSegregateAggregates)
      PFact->SetFactory("A", segregatedAFact);  // make sure that prolongator smoothing does not
                                                // disturb segregation of transfer operators
    // PFact->SetFactory("A",singleNodeAFact);
    // PFact->SetFactory("A",slaveTransferAFactory);  // produces nans
    RFact = Teuchos::rcp(new GenericRFactory());
  }

  // define nullspace factory AFTER tentative PFactory (that generates the nullspace for the coarser
  // levels) use same nullspace factory for all multigrid levels therefor we have to create one
  // instance of NullspaceFactory and use it for all FactoryManager objects (note: here, we have one
  // FactoryManager object per level)
  Teuchos::RCP<NullspaceFactory> nspFact = Teuchos::rcp(new NullspaceFactory("Nullspace"));
  nspFact->SetFactory("Nullspace", PtentFact);

  // RAP factory with inter-level transfer of segregation block information (map extractor)
  Teuchos::RCP<RAPFactory> AcFact = Teuchos::rcp(new RAPFactory());
  AcFact->SetFactory("P", PFact);
  AcFact->SetFactory("R", RFact);
  AcFact->SetParameter("RepairMainDiagonal", Teuchos::ParameterEntry(true));

  // write out aggregates
#if 0
  Teuchos::RCP<MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> > aggExpFact = Teuchos::rcp(new MueLu::AggregationExportFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>());
  aggExpFact->SetParameter("Output filename",Teuchos::ParameterEntry(std::string("aggs_%TIMESTEP(%ITER)_level%LEVELID_proc%PROCID.out")));
  if(params.isSublist("Linear System properties")) {
      const Teuchos::ParameterList & linSystemProps = params.sublist("Linear System properties");
      //epMasterDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map> > ("contact masterDofMap");
      aggExpFact->SetParameter("Output file: time step",Teuchos::ParameterEntry(linSystemProps.get< int > ("time step")));
      aggExpFact->SetParameter("Output file: iter",Teuchos::ParameterEntry(linSystemProps.get< int > ("iter")));
  }
  aggExpFact->SetFactory("Aggregates",UCAggFact);
  aggExpFact->SetFactory("DofsPerNode",dropFact);
  AcFact->AddTransferFactory(aggExpFact);
#endif

  // transfer maps to coarser grids
  if (xSingleNodeAggMap != Teuchos::null)
  {
    Teuchos::RCP<MapTransferFactory> cmTransFact = Teuchos::rcp(new MapTransferFactory());
    cmTransFact->SetParameter("map: factory", Teuchos::ParameterEntry(std::string("NoFactory")));
    cmTransFact->SetParameter(
        "map: name", Teuchos::ParameterEntry(std::string("SingleNodeAggDofMap")));
    cmTransFact->SetFactory("P", PtentFact);
    AcFact->AddTransferFactory(cmTransFact);
  }
  if (xSlaveDofMap != Teuchos::null)
  {  // needed for ContactAFilterFactory, IterationAFactory
    Teuchos::RCP<MapTransferFactory> cmTransFact = Teuchos::rcp(new MapTransferFactory());
    cmTransFact->SetParameter("map: factory", Teuchos::ParameterEntry(std::string("NoFactory")));
    cmTransFact->SetParameter("map: name", Teuchos::ParameterEntry(std::string("SlaveDofMap")));
    cmTransFact->SetFactory("P", PtentFact);
    AcFact->AddTransferFactory(cmTransFact);
  }
  if (xMasterDofMap != Teuchos::null)
  {  // needed for ContactAFilterFactory
    Teuchos::RCP<MapTransferFactory> cmTransFact = Teuchos::rcp(new MapTransferFactory());
    cmTransFact->SetParameter("map: factory", Teuchos::ParameterEntry(std::string("NoFactory")));
    cmTransFact->SetParameter("map: name", Teuchos::ParameterEntry(std::string("MasterDofMap")));
    cmTransFact->SetFactory("P", PtentFact);
    AcFact->AddTransferFactory(cmTransFact);
  }
  if (xNearZeroDiagMap != Teuchos::null)
  {
    Teuchos::RCP<MapTransferFactory> cmTransFact = Teuchos::rcp(new MapTransferFactory());
    cmTransFact->SetParameter("map: factory", Teuchos::ParameterEntry(std::string("NoFactory")));
    cmTransFact->SetParameter("map: name", Teuchos::ParameterEntry(std::string("NearZeroDiagMap")));
    cmTransFact->SetFactory("P", PtentFact);
    AcFact->AddTransferFactory(cmTransFact);
  }

  ///////////////////////////////////////////////////////////////////////
  // introduce rebalancing
  ///////////////////////////////////////////////////////////////////////
#ifdef HAVE_MUELU_ISORROPIA
  Teuchos::RCP<Factory> RebalancedPFact = Teuchos::null;
  Teuchos::RCP<Factory> RebalancedRFact = Teuchos::null;
  Teuchos::RCP<Factory> RepartitionFact = Teuchos::null;
  Teuchos::RCP<RebalanceAcFactory> RebalancedAFact = Teuchos::null;
  if (bDoRepartition)
  {
    // The Factory Manager will be configured to return the rebalanced versions of P, R, A by
    // default. Everytime we want to use the non-rebalanced versions, we need to explicitly define
    // the generating factory.
    RFact->SetFactory("P", PFact);
    //
    AcFact->SetFactory("P", PFact);
    AcFact->SetFactory("R", RFact);

    // define rebalancing factory for coarse matrix
    Teuchos::RCP<AmalgamationFactory> rebAmalgFact = Teuchos::rcp(new AmalgamationFactory());
    rebAmalgFact->SetFactory("A", AcFact);

    // create amalgamated "Partition"
    Teuchos::RCP<MueLu::IsorropiaInterface<LO, GO, NO>> isoInterface =
        Teuchos::rcp(new MueLu::IsorropiaInterface<LO, GO, NO>());
    isoInterface->SetFactory("A", AcFact);
    isoInterface->SetFactory("UnAmalgamationInfo", rebAmalgFact);

    // create "Partition" by unamalgamtion
    Teuchos::RCP<MueLu::RepartitionInterface<LO, GO, NO>> repInterface =
        Teuchos::rcp(new MueLu::RepartitionInterface<LO, GO, NO>());
    repInterface->SetFactory("A", AcFact);
    repInterface->SetFactory("AmalgamatedPartition", isoInterface);
    // repInterface->SetFactory("UnAmalgamationInfo", rebAmalgFact); // not necessary?

    // Repartitioning (creates "Importer" from "Partition")
    RepartitionFact = Teuchos::rcp(new RepartitionFactory());
    {
      Teuchos::ParameterList paramList;
      paramList.set("minRowsPerProcessor", /*10*/ optMinRowsPerProc);
      paramList.set("nonzeroImbalance", /*1.2*/ optNnzImbalance);
      RepartitionFact->SetParameterList(paramList);
    }
    RepartitionFact->SetFactory("A", AcFact);
    RepartitionFact->SetFactory("Partition", repInterface);

    if (reuseStrategy == "Ptent")
      h->Keep("Importer", RepartitionFact.get());  // keep Importer from RepartitionFact!!

    // Reordering of the transfer operators
    RebalancedPFact = Teuchos::rcp(new RebalanceTransferFactory());
    RebalancedPFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Interpolation")));
    RebalancedPFact->SetFactory("P", PFact);

    RebalancedRFact = Teuchos::rcp(new RebalanceTransferFactory());
    RebalancedRFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Restriction")));
    RebalancedRFact->SetFactory("R", RFact);
    RebalancedRFact->SetFactory("Nullspace", PtentFact);
    if (reuseStrategy == "Ptent") RebalancedRFact->DisableMultipleCheckGlobally();

    // Compute Ac from rebalanced P and R
    RebalancedAFact = Teuchos::rcp(new RebalanceAcFactory());
    RebalancedAFact->SetFactory("A", AcFact);

    // Rebalance maps
    Teuchos::RCP<RebalanceMapFactory> rebFact = Teuchos::rcp(new RebalanceMapFactory());
    rebFact->SetParameter("Map name", Teuchos::ParameterEntry(std::string("SingleNodeAggDofMap")));
    RebalancedAFact->AddRebalanceFactory(rebFact);

    Teuchos::RCP<RebalanceMapFactory> rebFact2 = Teuchos::rcp(new RebalanceMapFactory());
    rebFact2->SetParameter("Map name", Teuchos::ParameterEntry(std::string("SlaveDofMap")));
    RebalancedAFact->AddRebalanceFactory(rebFact2);

    Teuchos::RCP<RebalanceMapFactory> rebFact3 = Teuchos::rcp(new RebalanceMapFactory());
    rebFact3->SetParameter("Map name", Teuchos::ParameterEntry(std::string("MasterDofMap")));
    RebalancedAFact->AddRebalanceFactory(rebFact3);

    Teuchos::RCP<RebalanceMapFactory> rebFact4 = Teuchos::rcp(new RebalanceMapFactory());
    rebFact4->SetParameter("Map name", Teuchos::ParameterEntry(std::string("NearZeroDiagMap")));
    RebalancedAFact->AddRebalanceFactory(rebFact4);
  }
#endif  // #ifdef HAVE_MUELU_ISORROPIA

  ///////////////////////////////////////////////////////////////////////
  // prepare factory managers
  ///////////////////////////////////////////////////////////////////////


  // vecManager_.reserve(maxLevels);
  for (int i = 0; i < maxLevels; i++)
  {
    Teuchos::ParameterList pp(params);

    vecManager_.push_back(Teuchos::rcp(new FactoryManager()));

    vecManager_[i]->SetFactory("Aggregates", UCAggFact);
    vecManager_[i]->SetFactory("Graph", dropFact);
    vecManager_[i]->SetFactory("DofsPerNode", dropFact);
#ifdef HAVE_MUELU_ISORROPIA
    if (bDoRepartition)
    {
      vecManager_[i]->SetFactory("A", RebalancedAFact);
      vecManager_[i]->SetFactory("P", RebalancedPFact);
      vecManager_[i]->SetFactory("R", RebalancedRFact);
      vecManager_[i]->SetFactory("Nullspace", RebalancedRFact);
      vecManager_[i]->SetFactory("Importer", RepartitionFact);
      vecManager_[i]->SetFactory(
          "Ptent", PtentFact);  // same prolongator and restrictor factories for all levels
    }
    else
    {
#endif  // #ifdef HAVE_MUELU_ISORROPIA
      vecManager_[i]->SetFactory(
          "Nullspace", nspFact);  // use same nullspace factory throughout all multigrid levels
      vecManager_[i]->SetFactory("A", AcFact);  // same RAP factory for all levels
      vecManager_[i]->SetFactory(
          "P", PFact);  // same prolongator and restrictor factories for all levels
      vecManager_[i]->SetFactory(
          "Ptent", PtentFact);  // same prolongator and restrictor factories for all levels
      vecManager_[i]->SetFactory(
          "R", RFact);  // same prolongator and restrictor factories for all levels
#ifdef HAVE_MUELU_ISORROPIA
    }
#endif
  }

  return h;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>
LINALG::SOLVER::MueLuContactPreconditioner2::SetupSmoothers(
    const Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>& h, const Teuchos::ParameterList& params)
{
  // extract additional maps from parameter list
  // these maps are provided by the STR::TimInt::PrepareContactMeshtying routine, that
  // has access to the contact manager class
  Teuchos::RCP<Map> xSingleNodeAggMap = Teuchos::null;
  if (params.isSublist("Linear System properties"))
  {
    const Teuchos::ParameterList& linSystemProps = params.sublist("Linear System properties");
    if (linSystemProps.isParameter("non diagonal-dominant row map"))
      xSingleNodeAggMap = linSystemProps.get<Teuchos::RCP<Map>>("non diagonal-dominant row map");
  }

  // for the Jacobi/SGS smoother we wanna change the input matrix A and set Dirichlet bcs for the
  // (active?) slave dofs rebuild a new singleNodeAFact
  h->GetLevel(0)->Delete("A", singleNodeAFact_.get());
  singleNodeAFact_ = Teuchos::null;  // delete the old singleNodeAFact_ if there was one before
  if (xSingleNodeAggMap != Teuchos::null)
  {
    singleNodeAFact_ =
        Teuchos::rcp(new MueLu::IterationAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
            /*"SingleNodeAggDofMap",MueLu::NoFactory::getRCP()*/));
    singleNodeAFact_->SetParameter(
        "map: name", Teuchos::ParameterEntry(std::string("SingleNodeAggDofMap")));
    singleNodeAFact_->SetParameter(
        "map: factory", Teuchos::ParameterEntry(std::string("NoFactory")));
    // keep singleNodeAFact since it's needed in the solution phase by MyTrilinosSmoother
    h->GetLevel(0)->Keep("A", singleNodeAFact_.get());
  }

  // set Smoother and CoarseSolver factories in FactoryManager
  int maxLevels = 10;
  if (params.isParameter("max levels")) maxLevels = params.get<int>("max levels");

  // coarse level smoother/solver
  Teuchos::RCP<SmootherFactory> coarsestSmooFact;
  coarsestSmooFact = GetContactCoarsestSolverFactory(
      params, Teuchos::null);  // use full matrix A on coarsest level (direct solver)

  // vecManager_.reserve(maxLevels);
  for (int i = 0; i < maxLevels; i++)
  {
    Teuchos::ParameterList pp(params);

    // fine/intermedium level smoother
    Teuchos::RCP<SmootherFactory> SmooFactFine = Teuchos::null;
    if (xSingleNodeAggMap != Teuchos::null)
      SmooFactFine = GetContactSmootherFactory(
          pp, i, singleNodeAFact_);  // use filtered matrix on fine and intermedium levels
    else
      SmooFactFine = GetContactSmootherFactory(pp, i, Teuchos::null);

    // Configure FactoryManager
    if (SmooFactFine != Teuchos::null)
      vecManager_[i]->SetFactory("Smoother",
          SmooFactFine);  // Hierarchy.Setup uses TOPSmootherFactory, that only needs "Smoother"
    vecManager_[i]->SetFactory("CoarseSolver", coarsestSmooFact);
  }

  return h;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>
LINALG::SOLVER::MueLuContactPreconditioner2::SetupHierarchy(
    const Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>& h, const Teuchos::ParameterList& params,
    const Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>>& A,
    const Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nsp)
{
  bool bIsLastLevel = false;
  int maxLevels = 10;
  if (params.isParameter("max levels")) maxLevels = params.get<int>("max levels");

  ////////////////////////////

  // use new Hierarchy::Setup routine
  if (maxLevels == 1)
  {
    bIsLastLevel =
        h->Setup(0, Teuchos::null, vecManager_[0], Teuchos::null);  // 1 level "multigrid" method
  }
  else
  {
    bIsLastLevel =
        h->Setup(0, Teuchos::null, vecManager_[0], vecManager_[1]);  // first (finest) level
    for (int i = 1; i < maxLevels - 1; i++)
    {  // intermedium levels
      if (bIsLastLevel == true) break;
      bIsLastLevel = h->Setup(i, vecManager_[i - 1], vecManager_[i], vecManager_[i + 1]);
    }
    if (bIsLastLevel == false)
    {  // coarsest level
      bIsLastLevel = h->Setup(
          maxLevels - 1, vecManager_[maxLevels - 2], vecManager_[maxLevels - 1], Teuchos::null);
    }
  }

  return h;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactPreconditioner2::GetContactSmootherFactory(
    const Teuchos::ParameterList& paramList, int level, const Teuchos::RCP<FactoryBase>& AFact)
{
  char levelchar[11];
  sprintf(levelchar, "(level %d)", level);
  std::string levelstr(levelchar);

  if (paramList.isSublist("smoother: list " + levelstr) == false) return Teuchos::null;
  TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist("smoother: list " + levelstr) == false,
      MueLu::Exceptions::RuntimeError,
      "MueLu::Interpreter: no ML smoother parameter list for level. error.");

  std::string type =
      paramList.sublist("smoother: list " + levelstr).get<std::string>("smoother: type");
  TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), MueLu::Exceptions::RuntimeError,
      "MueLu::Interpreter: no ML smoother type for level. error.");

  const Teuchos::ParameterList smolevelsublist = paramList.sublist("smoother: list " + levelstr);

  Teuchos::RCP<MueLu::SmootherPrototype<SC, LO, GO, NO>> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  Teuchos::RCP<MueLu::SmootherFactory<SC, LO, GO, NO>> SmooFact;

  if (type == "Jacobi")
  {
    if (smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", smolevelsublist.get<int>("smoother: sweeps"));
    if (smolevelsublist.get<double>("smoother: damping factor"))
      ifpackList.set(
          "relaxation: damping factor", smolevelsublist.get<double>("smoother: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Jacobi");
    smooProto =
        Teuchos::rcp(new MueLu::MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
            "SingleNodeAggDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact));
  }
  else if (type == "Gauss-Seidel")
  {
    if (smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", smolevelsublist.get<int>("smoother: sweeps"));
    if (smolevelsublist.get<double>("smoother: damping factor"))
      ifpackList.set(
          "relaxation: damping factor", smolevelsublist.get<double>("smoother: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Gauss-Seidel");
    smooProto =
        Teuchos::rcp(new MueLu::MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
            "SingleNodeAggDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact));
  }
  else if (type == "symmetric Gauss-Seidel")
  {
    if (smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", smolevelsublist.get<int>("smoother: sweeps"));
    if (smolevelsublist.get<double>("smoother: damping factor"))
      ifpackList.set(
          "relaxation: damping factor", smolevelsublist.get<double>("smoother: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto =
        Teuchos::rcp(new MueLu::MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
            "SingleNodeAggDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact));
    // std::cout << "built symm GS: " << smooProto << std::endl;
  }
  else if (type == "Chebyshev")
  {
    ifpackType = "CHEBYSHEV";
    if (smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set("chebyshev: degree", smolevelsublist.get<int>("smoother: sweeps"));
    smooProto =
        Teuchos::rcp(new MueLu::MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
            "SingleNodeAggDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, 0, AFact));
    // TODO what about the other parameters
  }
  else if (type == "IFPACK")
  {
#ifdef HAVE_MUELU_IFPACK
    // TODO change to TrilinosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
    ifpackType =
        paramList.sublist("smoother: list " + levelstr).get<std::string>("smoother: ifpack type");
    if (ifpackType == "ILU")
    {
      ifpackList.set<int>("fact: level-of-fill",
          (int)smolevelsublist.get<double>("smoother: ifpack level-of-fill"));
      ifpackList.set("partitioner: overlap", smolevelsublist.get<int>("smoother: ifpack overlap"));
      int overlap = smolevelsublist.get<int>("smoother: ifpack overlap");
      // smooProto = MueLu::GetIfpackSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node>(ifpackType,
      // ifpackList,smolevelsublist.get<int>("smoother: ifpack overlap"),AFact);
      smooProto =
          Teuchos::rcp(new MueLu::MyTrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
              "SingleNodeAggDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, overlap,
              AFact));
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
          "MueLu::Interpreter: unknown ML smoother type " + type +
              " (IFPACK) not supported by MueLu. Only ILU is supported.");
#else   // HAVE_MUELU_IFPACK
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
        "MueLu::Interpreter: MueLu compiled without Ifpack support");
#endif  // HAVE_MUELU_IFPACK
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
        "MueLu::Interpreter: unknown ML smoother type " + type + " not supported by MueLu.");
  }

  // create smoother factory
  SmooFact = Teuchos::rcp(new MueLu::SmootherFactory<SC, LO, GO, NO>(smooProto));

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if (smolevelsublist.isParameter("smoother: pre or post"))
    preorpost = smolevelsublist.get<std::string>("smoother: pre or post");

  if (preorpost == "pre")
  {
    SmooFact->SetSmootherPrototypes(smooProto, Teuchos::null);
  }
  else if (preorpost == "post")
  {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smooProto);
  }

  return SmooFact;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactPreconditioner2::GetContactCoarsestSolverFactory(
    const Teuchos::ParameterList& paramList, const Teuchos::RCP<FactoryBase>& AFact)
{
  std::string type = "";  // use default defined by AmesosSmoother or Amesos2Smoother

  if (paramList.isParameter("coarse: type")) type = paramList.get<std::string>("coarse: type");

  Teuchos::RCP<MueLu::SmootherPrototype<SC, LO, GO, NO>> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  Teuchos::RCP<MueLu::SmootherFactory<SC, LO, GO, NO>> SmooFact;

  if (type == "Jacobi")
  {
    if (paramList.isParameter("coarse: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", paramList.get<int>("coarse: sweeps"));
    else
      ifpackList.set<int>("relaxation: sweeps", 1);
    if (paramList.isParameter("coarse: damping factor"))
      ifpackList.set("relaxation: damping factor", paramList.get<double>("coarse: damping factor"));
    else
      ifpackList.set("relaxation: damping factor", 1.0);
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Jacobi");
    smooProto = rcp(new MueLu::TrilinosSmoother<SC, LO, GO, NO>(ifpackType, ifpackList, 0));
    smooProto->SetFactory("A", AFact);
  }
  else if (type == "Gauss-Seidel")
  {
    if (paramList.isParameter("coarse: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", paramList.get<int>("coarse: sweeps"));
    else
      ifpackList.set<int>("relaxation: sweeps", 1);
    if (paramList.isParameter("coarse: damping factor"))
      ifpackList.set("relaxation: damping factor", paramList.get<double>("coarse: damping factor"));
    else
      ifpackList.set("relaxation: damping factor", 1.0);
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Gauss-Seidel");
    smooProto = rcp(new MueLu::TrilinosSmoother<SC, LO, GO, NO>(ifpackType, ifpackList, 0));
    smooProto->SetFactory("A", AFact);
  }
  else if (type == "symmetric Gauss-Seidel")
  {
    if (paramList.isParameter("coarse: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", paramList.get<int>("coarse: sweeps"));
    else
      ifpackList.set<int>("relaxation: sweeps", 1);
    if (paramList.isParameter("coarse: damping factor"))
      ifpackList.set("relaxation: damping factor", paramList.get<double>("coarse: damping factor"));
    else
      ifpackList.set("relaxation: damping factor", 1.0);
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto = rcp(new MueLu::TrilinosSmoother<SC, LO, GO, NO>(ifpackType, ifpackList, 0));
    smooProto->SetFactory("A", AFact);
  }
  else if (type == "Chebyshev")
  {
    ifpackType = "CHEBYSHEV";
    if (paramList.isParameter("coarse: sweeps"))
      ifpackList.set("chebyshev: degree", paramList.get<int>("coarse: sweeps"));
    if (paramList.isParameter("coarse: Chebyshev alpha"))
      ifpackList.set("chebyshev: alpha", paramList.get<double>("coarse: Chebyshev alpha"));
    smooProto = rcp(new MueLu::TrilinosSmoother<SC, LO, GO, NO>(ifpackType, ifpackList, 0));
    smooProto->SetFactory("A", AFact);
  }
  else if (type == "IFPACK")
  {
#ifdef HAVE_MUELU_IFPACK
    // TODO change to TrilinosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
    /* ifpackType = paramList.get<std::string>("coarse: ifpack type");
     if(ifpackType == "ILU") {
       ifpackList.set<int>("fact: level-of-fill", (int)paramList.get<double>("coarse: ifpack
     level-of-fill")); ifpackList.set("partitioner: overlap", paramList.get<int>("coarse: ifpack
     overlap")); smooProto =
     MueLu::GetIfpackSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node>(ifpackType, ifpackList,
     paramList.get<int>("coarse: ifpack overlap"), AFact);
     }
     else*/
    //  TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::Interpreter:
    //  unknown ML smoother type " + type + " (IFPACK) not supported by MueLu. Only ILU is
    //  supported.");
    ifpackType = paramList.get<std::string>("coarse: ifpack type");
    if (ifpackType == "ILU")
    {
      ifpackList.set<int>(
          "fact: level-of-fill", (int)paramList.get<double>("coarse: ifpack level-of-fill"));
      ifpackList.set("partitioner: overlap", paramList.get<int>("coarse: ifpack overlap"));
      // int overlap = paramList.get<int>("coarse: ifpack overlap");
      // smooProto = Teuchos::rcp( new MueLu::PermutingSmoother<Scalar,LocalOrdinal, GlobalOrdinal,
      // Node>("SingleNodeAggDofMap", MueLu::NoFactory::getRCP(), ifpackType, ifpackList, overlap,
      // AFact) );
      smooProto = MueLu::GetIfpackSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
          ifpackType, ifpackList, paramList.get<int>("coarse: ifpack overlap"));
      smooProto->SetFactory("A", AFact);
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
          "MueLu::Interpreter: unknown ML smoother type " + type +
              " (IFPACK) not supported by MueLu. Only ILU is supported.");

#else   // HAVE_MUELU_IFPACK
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
        "MueLu::Interpreter: MueLu compiled without Ifpack support");
#endif  // HAVE_MUELU_IFPACK
  }
  else if (type == "Amesos-Superlu")
  {
    smooProto =
        Teuchos::rcp(new MueLu::DirectSolver<SC, LO, GO, NO>("Superlu", Teuchos::ParameterList()));
    smooProto->SetFactory("A", AFact);
  }
  else if (type == "Amesos-Superludist")
  {
    smooProto = Teuchos::rcp(
        new MueLu::DirectSolver<SC, LO, GO, NO>("Superludist", Teuchos::ParameterList()));
    smooProto->SetFactory("A", AFact);
  }
  else if (type == "Amesos-KLU")
  {
    smooProto =
        Teuchos::rcp(new MueLu::DirectSolver<SC, LO, GO, NO>("Klu", Teuchos::ParameterList()));
    smooProto->SetFactory("A", AFact);
  }
  else if (type == "Amesos-UMFPACK")
  {
    smooProto =
        Teuchos::rcp(new MueLu::DirectSolver<SC, LO, GO, NO>("Umfpack", Teuchos::ParameterList()));
    smooProto->SetFactory("A", AFact);
  }
  else if (type == "")
  {
    smooProto = Teuchos::rcp(new MueLu::DirectSolver<SC, LO, GO, NO>("", Teuchos::ParameterList()));
    smooProto->SetFactory("A", AFact);
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
        "MueLu::Interpreter: unknown coarsest solver type. '" << type
                                                              << "' not supported by MueLu.");
  }

  // create smoother factory
  TEUCHOS_TEST_FOR_EXCEPTION(smooProto == Teuchos::null, MueLu::Exceptions::RuntimeError,
      "MueLu::Interpreter: no smoother prototype. fatal error.");
  SmooFact = rcp(new MueLu::SmootherFactory<SC, LO, GO, NO>(smooProto));

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if (paramList.isParameter("coarse: pre or post"))
    preorpost = paramList.get<std::string>("coarse: pre or post");

  if (preorpost == "pre")
  {
    SmooFact->SetSmootherPrototypes(smooProto, Teuchos::null);
  }
  else if (preorpost == "post")
  {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smooProto);
  }

  return SmooFact;
}

#endif
