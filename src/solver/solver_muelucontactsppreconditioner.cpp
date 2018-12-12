/*!----------------------------------------------------------------------
\file solver_muelucontactsppreconditioner.cpp

\brief Implementation of MueLu-base contact saddle-point preconditioners

\level 2

\maintainer Matthias Mayr
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
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_CrsMatrix.hpp>  // for merging blocked operator

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
#include <MueLu_BlockedRAPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_NullspaceFactory.hpp>

#include <MueLu_Aggregates.hpp>
#include <MueLu_AggregationAlgorithmBase.hpp>
#include <MueLu_AggregationExportFactory.hpp>
#include <MueLu_BlockedPFactory.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_SchurComplementFactory.hpp>
#include <MueLu_BraessSarazinSmoother.hpp>
#include <MueLu_SimpleSmoother.hpp>
#include <MueLu_IndefBlockedDiagonalSmoother.hpp>
#include <MueLu_UzawaSmoother.hpp>
#include <MueLu_PermutingSmoother.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_MapTransferFactory.hpp>
#include <MueLu_BlockedCoarseMapFactory.hpp>
#include <MueLu_PermutationFactory.hpp>

#include <MueLu_MLParameterListInterpreter.hpp>

#ifdef HAVE_MUELU_ISORROPIA
#include "MueLu_IsorropiaInterface.hpp"
#include "MueLu_RepartitionInterface.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RebalanceBlockInterpolationFactory.hpp"
#include "MueLu_RebalanceBlockRestrictionFactory.hpp"
#include "MueLu_RebalanceBlockAcFactory.hpp"
#include "MueLu_RebalanceMapFactory.hpp"
#endif

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include <MueLu_EpetraOperator.hpp>  // Aztec interface

#include "muelu/MueLu_ContactTransferFactory_decl.hpp"
#include "muelu/MueLu_ContactASlaveDofFilterFactory_decl.hpp"
#include "muelu/MueLu_ContactSPAggregationFactory_decl.hpp"
#include "muelu/MueLu_MyTrilinosSmoother_decl.hpp"
#include "muelu/MueLu_MeshtyingSPAmalgamationFactory_decl.hpp"
#include "muelu/MueLu_SelectiveSaPFactory_decl.hpp"
#include "muelu/MueLu_ContactSPRepartitionInterface_decl.hpp"

#include "solver_muelucontactsppreconditioner.H"

typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixWrap;
typedef Xpetra::EpetraCrsMatrix EpetraCrsMatrix;
typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
typedef Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> StridedMap;
typedef Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node> StridedMapFactory;
typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedCrsMatrix;

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
typedef MueLu::SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> SubBlockAFactory;
typedef MueLu::PFactory PFactory;
typedef MueLu::TwoLevelFactoryBase TwoLevelFactoryBase;
typedef MueLu::CoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CoarseMapFactory;
typedef MueLu::BlockedCoarseMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    BlockedCoarseMapFactory;
typedef MueLu::RebalanceBlockRestrictionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    RebalanceBlockRestrictionFactory;
typedef MueLu::UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node> UzawaSmoother;
typedef MueLu::IndefBlockedDiagonalSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    IndefBlockedDiagonalSmoother;
typedef MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> SmootherPrototype;
typedef MueLu::PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> PermutationFactory;
typedef MueLu::SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    SchurComplementFactory;
typedef MueLu::DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node> DirectSolver;
typedef MueLu::BlockedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedPFactory;
typedef MueLu::BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedRAPFactory;
typedef MueLu::RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    RebalanceBlockAcFactory;
typedef MueLu::SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node> SimpleSmoother;
typedef MueLu::BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    BraessSarazinSmoother;
typedef MueLu::TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node> TrilinosSmoother;
typedef MueLu::RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    RebalanceBlockInterpolationFactory;
typedef MueLu::PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node> PermutingSmoother;
typedef MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node> Utils;

typedef Scalar SC;
typedef LocalOrdinal LO;
typedef GlobalOrdinal GO;
typedef Node NO;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuContactSpPreconditioner::MueLuContactSpPreconditioner(
    FILE* outfile, Teuchos::ParameterList& mllist)
    : PreconditionerType(outfile), mllist_(mllist)
{
  P_ = Teuchos::null;
  Pmatrix_ = Teuchos::null;
  H_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuContactSpPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  // Check whether input matrix is an actual blocked operator
  Teuchos::RCP<BlockSparseMatrixBase> A =
      Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) dserror("matrix is not a BlockSparseMatrix");

  // store blocked operator
  Pmatrix_ = A;

  ///////////////////////////////////////////////////////////////////////
  // Some preparations common for re-creating and re-using the preconditioner
  ///////////////////////////////////////////////////////////////////////

  // prepare nullspace vector for MueLu (block A11 only)
  const int numdf = mllist_.get<int>("PDE equations", -1);
  const int dimns = mllist_.get<int>("null space: dimension", -1);
  if (numdf == -1) dserror("Multigrid parameter 'PDE equations' wrong. It has to be > 0.");
  if (dimns == -1) dserror("Multigrid parameter 'null space: dimension' wrong. It has to be > 0.");

  // create a Teuchos::Comm from EpetraComm
  Teuchos::RCP<const Teuchos::Comm<int>> comm = Xpetra::toXpetra(A->RangeMap(0).Comm());

  /* Extract additional maps from parameter list
   *
   * These maps are provided by the STR::TimInt::PrepareContactMeshtying routine, that has access
   * to the contact manager class
   *
   * Note: Baci provides Epetra_Map objects. We will transform them to Xpetra::Map later.
   */
  Teuchos::RCP<Epetra_Map> epMasterDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epSlaveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epActiveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epInnerDofMap = Teuchos::null;
  Teuchos::RCP<Map> xSingleNodeAggMap = Teuchos::null;
  Teuchos::RCP<Map> xNearZeroDiagMap = Teuchos::null;
  if (mllist_.isSublist("Linear System properties"))
  {
    const Teuchos::ParameterList& linSystemProps = mllist_.sublist("Linear System properties");
    // extract information provided by solver (e.g. PermutedAztecSolver)
    epMasterDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact masterDofMap");
    epSlaveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact slaveDofMap");
    epActiveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact activeDofMap");
    epInnerDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact innerDofMap");
  }

  ///////////////////////////////////////////////////////////////////////
  // Transform from epetra to xpetra
  ///////////////////////////////////////////////////////////////////////

  // Transform maps
  Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap = Teuchos::rcp(new Xpetra::EpetraMap(epSlaveDofMap));

  // Get maps and matrix blocks for blocked operator
  Teuchos::RCP<const Map> fullrangemap = Teuchos::null;
  Teuchos::RCP<CrsMatrix> xA11 = Teuchos::null;
  Teuchos::RCP<CrsMatrix> xA12 = Teuchos::null;
  Teuchos::RCP<CrsMatrix> xA21 = Teuchos::null;
  Teuchos::RCP<CrsMatrix> xA22 = Teuchos::null;
  if (create)
  {
    // (Re-)create the preconditioner, so extract everything from the input matrix 'A'.

    // Prepare maps for blocked operator
    fullrangemap = Teuchos::rcp(new Xpetra::EpetraMap(Teuchos::rcpFromRef(A->FullRangeMap())));

    // Transform matrix blocks
    xA11 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 0).EpetraMatrix()));
    xA12 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 1).EpetraMatrix()));
    xA21 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 0).EpetraMatrix()));
    xA22 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 1).EpetraMatrix()));
  }
  else
  {
    // Re-use the preconditioner, so extract everything from the existing preconditioner 'Pmatrix_'.

    // create maps
    fullrangemap =
        Teuchos::rcp(new Xpetra::EpetraMap(Teuchos::rcpFromRef(Pmatrix_->FullRangeMap())));

    xA11 = Teuchos::rcp(new EpetraCrsMatrix(Pmatrix_->Matrix(0, 0).EpetraMatrix()));
    xA12 = Teuchos::rcp(new EpetraCrsMatrix(Pmatrix_->Matrix(0, 1).EpetraMatrix()));
    xA21 = Teuchos::rcp(new EpetraCrsMatrix(Pmatrix_->Matrix(1, 0).EpetraMatrix()));
    xA22 = Teuchos::rcp(new EpetraCrsMatrix(Pmatrix_->Matrix(1, 1).EpetraMatrix()));
  }

  /* Define strided maps
   *
   * We have 'numdf' Lagrange multipliers per node at the contact interface,
   * so we also add 'numdf' entries to map of the Lagrange multipliers.
   */
  std::vector<size_t> stridingInfo1;
  std::vector<size_t> stridingInfo2;
  stridingInfo1.push_back(numdf);
  stridingInfo2.push_back(numdf);
  Teuchos::RCP<StridedMap> strMap1 = Teuchos::rcp(new StridedMap(xA11->getRowMap(), stridingInfo1,
      xA11->getRowMap()->getIndexBase(), -1 /* stridedBlock */, 0 /*globalOffset*/));
  Teuchos::RCP<StridedMap> strMap2 = Teuchos::rcp(new StridedMap(xA22->getRowMap(), stridingInfo2,
      xA22->getRowMap()->getIndexBase(), -1 /* stridedBlock */, 0 /*globalOffset*/));

  // build map extractor
  std::vector<Teuchos::RCP<const Map>> xmaps;
  xmaps.push_back(strMap1);
  xmaps.push_back(strMap2);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO>> map_extractor =
      Xpetra::MapExtractorFactory<Scalar, LO, GO>::Build(fullrangemap, xmaps);

  // build blocked Xpetra operator
  Teuchos::RCP<BlockedCrsMatrix> bOp =
      Teuchos::rcp(new BlockedCrsMatrix(map_extractor, map_extractor, 10));
  bOp->setMatrix(0, 0, xA11);
  bOp->setMatrix(0, 1, xA12);
  bOp->setMatrix(1, 0, xA21);
  bOp->setMatrix(1, 1, xA22);
  bOp->fillComplete();

  // Re-create or re-use the preconditioner?
  if (create)
  {
    // (Re-)create the preconditioner

    // free old matrix first
    P_ = Teuchos::null;

    // interpret ML parameters
    int maxLevels = 3;
    int verbosityLevel = 10;
    int maxCoarseSize = 100;
    int minPerAgg = 9;   // optimal for 2d
    int maxPerAgg = 27;  // optimal for 3d
    int maxNbrAlreadySelected = 0;
    std::string agg_type = "Uncoupled";
    double agg_damping = 0.0;
    if (mllist_.isParameter("max levels")) maxLevels = mllist_.get<int>("max levels");
    if (mllist_.isParameter("ML output")) verbosityLevel = mllist_.get<int>("ML output");
    if (mllist_.isParameter("coarse: max size"))
      maxCoarseSize = mllist_.get<int>("coarse: max size");
    if (mllist_.isParameter("aggregation: type"))
      agg_type = mllist_.get<std::string>("aggregation: type");
    if (mllist_.isParameter("aggregation: nodes per aggregate"))
      maxPerAgg = mllist_.get<int>("aggregation: nodes per aggregate");
    if (mllist_.isParameter("aggregation: min nodes per aggregate"))
      minPerAgg = mllist_.get<int>("aggregation: min nodes per aggregate");
    if (mllist_.isParameter("aggregation: damping factor"))
      agg_damping = mllist_.get<double>("aggregation: damping factor");

#ifdef HAVE_MUELU_ISORROPIA
    // parameters for repartitioning
    bool bDoRepartition = false;
    double optNnzImbalance = 1.3;
    int optMinRowsPerProc = 3000;

    if (mllist_.isParameter("muelu repartition: enable"))
    {
      if (mllist_.get<int>("muelu repartition: enable") == 1) bDoRepartition = true;
    }
    if (mllist_.isParameter("muelu repartition: max min ratio"))
      optNnzImbalance = mllist_.get<double>("muelu repartition: max min ratio");
    if (mllist_.isParameter("muelu repartition: min per proc"))
      optMinRowsPerProc = mllist_.get<int>("muelu repartition: min per proc");
#else
    dserror("Isorropia has not been compiled with Trilinos. Repartitioning is not working.");
#endif

    // translate verbosity parameter
    Teuchos::EVerbosityLevel eVerbLevel = Teuchos::VERB_NONE;
    if (verbosityLevel == 0) eVerbLevel = Teuchos::VERB_NONE;
    if (verbosityLevel > 0) eVerbLevel = Teuchos::VERB_LOW;
    if (verbosityLevel > 4) eVerbLevel = Teuchos::VERB_MEDIUM;
    if (verbosityLevel > 7) eVerbLevel = Teuchos::VERB_HIGH;
    if (verbosityLevel > 9) eVerbLevel = Teuchos::VERB_EXTREME;

    bool bIsMeshtying = false;  // assume a contact problem
    if (mllist_.isSublist("Linear System properties"))
    {
      const Teuchos::ParameterList& linSystemProps = mllist_.sublist("Linear System properties");

      if (linSystemProps.isParameter("ProblemType") &&
          linSystemProps.get<std::string>("ProblemType") == "meshtying")
        bIsMeshtying = true;
    }

    if (bIsMeshtying)
      dserror(
          "This is a contact preconditioner. Should we really end up here? It might be ok, but "
          "I'm "
          "not sure as we currently can't do saddle-point preconditioning for meshtying "
          "problems.");

    ///////////////////////////////////////////////////////////////////////
    // prepare nullspace for first block
    ///////////////////////////////////////////////////////////////////////

    // extract nullspace information from ML list
    Teuchos::RCP<MultiVector> nspVector11 =
        Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
            xA11->getRowMap(), dimns, true);
    Teuchos::RCP<std::vector<double>> nsdata =
        mllist_.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
    for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
    {
      Teuchos::ArrayRCP<Scalar> nspVector11i = nspVector11->getDataNonConst(i);
      const size_t myLength = nspVector11->getLocalLength();
      for (size_t j = 0; j < myLength; j++)
      {
        nspVector11i[j] = (*nsdata)[i * myLength + j];
      }
    }

    ///////////////////////////////////////////////////////////////////////
    // special aggregation strategy
    ///////////////////////////////////////////////////////////////////////

    // number of node rows (only displacement dofs)
    const LocalOrdinal nDofRows = strMap1->getNodeNumElements();

    // prepare aggCoarseStat
    // TODO rebuild node-based map
    // still problematic for repartitioning
    Teuchos::ArrayRCP<unsigned int> aggStat;
    if (nDofRows > 0) aggStat = Teuchos::arcp<unsigned int>(nDofRows / numdf);
    for (LocalOrdinal i = 0; i < nDofRows; ++i) aggStat[i / numdf] = MueLu::NodeState::READY;

    ///////////////////////////////////////////////////////////////////////
    // create Hierarchy
    ///////////////////////////////////////////////////////////////////////

    Teuchos::RCP<Hierarchy> H = rcp(new Hierarchy());
    H->setlib(Xpetra::UseEpetra);  // not very nice, but safe.
    H->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(eVerbLevel));
    H->SetMaxCoarseSize(Teuchos::as<Xpetra::global_size_t>(maxCoarseSize));
    H->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));
    H->GetLevel(0)->Set("Nullspace1", nspVector11);
    H->GetLevel(0)->Set("coarseAggStat", aggStat);
    H->GetLevel(0)->Set("SlaveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, Node>>(
                                           xSlaveDofMap));  // set map with active dofs

    Teuchos::RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
    Teuchos::RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
    A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
    A11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
    A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
    A22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

    ///////////////////////////////////////////////////////////////////////
    // set up block 11
    ///////////////////////////////////////////////////////////////////////

    Teuchos::RCP<AmalgamationFactory> amalgFact11 = Teuchos::rcp(new AmalgamationFactory());
    Teuchos::RCP<CoalesceDropFactory> dropFact11 = Teuchos::rcp(new CoalesceDropFactory());
    Teuchos::RCP<UncoupledAggregationFactory> UCAggFact11 =
        Teuchos::rcp(new UncoupledAggregationFactory());
    UCAggFact11->SetParameter(
        "aggregation: max selected neighbors", Teuchos::ParameterEntry(maxNbrAlreadySelected));
    UCAggFact11->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(minPerAgg));
    UCAggFact11->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(maxPerAgg));
    UCAggFact11->SetParameter(
        "aggregation: ordering", Teuchos::ParameterEntry(std::string("graph")));
    UCAggFact11->SetParameter(
        "aggregation: allow user-specified singletons", Teuchos::ParameterEntry(false));


    Teuchos::RCP<TentativePFactory> Ptent11Fact = Teuchos::rcp(new TentativePFactory());
    Teuchos::RCP<PFactory> P11Fact;
    Teuchos::RCP<TwoLevelFactoryBase> R11Fact;

    // Check for plain or smoothed aggregation
    if (agg_damping == 0.0)
    {
      // Plain aggregation: use tentative prolongation operator (PA-AMG)
      P11Fact = Ptent11Fact;
      R11Fact = Teuchos::rcp(new TransPFactory());
    }
    else if (agg_damping > 0.0)
    {
      // smoothed aggregation (SA-AMG)
      P11Fact =
          Teuchos::rcp(new MueLu::SelectiveSaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
      P11Fact->SetParameter("Damping factor", Teuchos::ParameterEntry(agg_damping));
      P11Fact->SetParameter("Damping strategy", Teuchos::ParameterEntry(std::string("User")));

      // feed SelectiveSAPFactory with information
      P11Fact->SetFactory("P", Ptent11Fact);
      P11Fact->SetFactory("A", A11Fact);
      /*
      // use user-given damping parameter
      P11Fact->SetParameter("Damping factor", Teuchos::ParameterEntry(agg_damping));
      P11Fact->SetParameter("Damping strategy", Teuchos::ParameterEntry(std::string("User")));
      // only prolongator smoothing for transfer operator basis functions which
      // correspond to non-slave rows in (permuted) matrix A
      // We use the tentative prolongator to detect the corresponding prolongator basis functions
      for given row gids.
      // Note: this ignores the permutations in A. In case, the matrix A has been permuted it can
      happen
      //       that problematic columns in Ptent are not corresponding to columns that belong to the
      //       with nonzero entries in slave rows. // TODO think more about this -> aggregation
      P11Fact->SetParameter("NonSmoothRowMapName",
      Teuchos::ParameterEntry(std::string("SlaveDofMap")));
      //PFact->SetParameter("NonSmoothRowMapName", Teuchos::ParameterEntry(std::string("")));
      P11Fact->SetFactory("NonSmoothRowMapFactory", MueLu::NoFactory::getRCP());

      // provide diagnostics of diagonal entries of current matrix A
      // if the solver object detects some significantly small entries on diagonal the contact
      // preconditioner can decide to skip transfer operator smoothing to increase robustness
      P11Fact->SetParameter("NearZeroDiagMapName",
      Teuchos::ParameterEntry(std::string("NearZeroDiagMap")));
      //PFact->SetParameter("NearZeroDiagMapName", Teuchos::ParameterEntry(std::string("")));
      P11Fact->SetFactory("NearZeroDiagMapFactory", MueLu::NoFactory::getRCP());
      */
      R11Fact = Teuchos::rcp(new GenericRFactory());
    }
    // ToDo (mayr.mt) Where do we set this? Is it clear that '-2.0' means PG-AMG with R=P^T?
    else if (agg_damping == -2.0)
    {
      // Petrov Galerkin PG-AMG smoothed aggregation (energy minimization in ML)
      P11Fact = Teuchos::rcp(new PgPFactory());
      P11Fact->SetFactory("P", Ptent11Fact);
      // PFact->SetFactory("A",singleNodeAFact);
      // PFact->SetFactory("A",slaveTransferAFactory);  // produces nans
      R11Fact = Teuchos::rcp(new TransPFactory());
      R11Fact->SetFactory("P", Ptent11Fact);
    }
    else
    {
      // Petrov Galerkin PG-AMG smoothed aggregation (energy minimization in ML)
      P11Fact = Teuchos::rcp(new PgPFactory());
      P11Fact->SetFactory("P", Ptent11Fact);
      // PFact->SetFactory("A",singleNodeAFact);
      // PFact->SetFactory("A",slaveTransferAFactory);  // produces nans
      R11Fact = Teuchos::rcp(new GenericRFactory());
    }

    Teuchos::RCP<NullspaceFactory> nspFact11 = Teuchos::rcp(new NullspaceFactory("Nullspace1"));
    nspFact11->SetFactory("Nullspace1", Ptent11Fact /*P11Fact*/);
    Teuchos::RCP<CoarseMapFactory> coarseMapFact11 = Teuchos::rcp(new CoarseMapFactory());
    coarseMapFact11->setStridingData(stridingInfo1);

    ///////////////////////////////////////////////////////////////////////
    // define factory manager for (1,1) block
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<FactoryManager> M11 = Teuchos::rcp(new FactoryManager());
    M11->SetFactory("A", A11Fact);
    M11->SetFactory("P", P11Fact);
    M11->SetFactory("Ptent", Ptent11Fact);
    M11->SetFactory("R", R11Fact);
    M11->SetFactory("Aggregates", UCAggFact11);
    M11->SetFactory("Nullspace", nspFact11);
    M11->SetFactory("Ptent", P11Fact);
    M11->SetFactory("Graph", dropFact11);
    M11->SetFactory("DofsPerNode", dropFact11);
    M11->SetFactory("CoarseMap", coarseMapFact11);
    M11->SetFactory("UnAmalgamationInfo", amalgFact11);
    M11->SetIgnoreUserData(true);  // always use data from factories defined in factory manager

    ///////////////////////////////////////////////////////////////////////
    // create default nullspace for block 2
    ///////////////////////////////////////////////////////////////////////

    const int dimNS2 = numdf;
    Teuchos::RCP<MultiVector> nspVector22 = MultiVectorFactory::Build(xA22->getRowMap(), dimNS2);

    for (int i = 0; i < dimNS2; ++i)
    {
      Teuchos::ArrayRCP<Scalar> nsValues22 = nspVector22->getDataNonConst(i);
      int numBlocks = nsValues22.size() / dimNS2;
      for (int j = 0; j < numBlocks; ++j)
      {
        nsValues22[j * dimNS2 + i] = 1.0;
      }
    }

    // set nullspace for block 2
    H->GetLevel(0)->Set("Nullspace2", nspVector22);

    ///////////////////////////////////////////////////////////////////////
    // set up block 2 factories
    ///////////////////////////////////////////////////////////////////////

    /* Use special aggregation routine which reconstructs aggregates for the Lagrange multipliers
     * using the aggregates for the displacements at the contact/meshtying interface to keep
     * correlation of interface nodes (displacements) and corresponding Lagrange multiplier DOFs
     */
    Teuchos::RCP<MueLu::ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
        UCAggFact22 = Teuchos::rcp(
            new MueLu::ContactSPAggregationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
                UCAggFact11, amalgFact11));

    // create amalgamation factory for (1,1) block depending on problem type (contact/meshtying)
    Teuchos::RCP<Factory> amalgFact22 = Teuchos::null;
    if (!bIsMeshtying)  // contact problem with (1,1) a non zero block
      amalgFact22 = Teuchos::rcp(new AmalgamationFactory());
    else  // meshtying problem with zero (1,1) block
      amalgFact22 = Teuchos::rcp(
          new MueLu::MeshtyingSPAmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
    amalgFact22->SetFactory("A", A22Fact);  // just to be sure, seems not to use factory from M22?

    // use tentative prolongation operator (prolongation operator smoothing doesn't make sense since
    // the A11 block is not valid for smoothing)
    Teuchos::RCP<TentativePFactory> P22Fact = Teuchos::rcp(new TentativePFactory());
    Teuchos::RCP<TransPFactory> R22Fact = Teuchos::rcp(new TransPFactory());
    R22Fact->SetFactory("P", P22Fact);
    Teuchos::RCP<NullspaceFactory> nspFact22 = Teuchos::rcp(new NullspaceFactory("Nullspace2"));
    nspFact22->SetFactory("Nullspace2", P22Fact);
    Teuchos::RCP<BlockedCoarseMapFactory> coarseMapFact22 =
        Teuchos::rcp(new BlockedCoarseMapFactory());
    coarseMapFact22->SetFactory(
        "CoarseMap", coarseMapFact11);  // feed BlockedCoarseMap with necessary input
    coarseMapFact22->SetFactory("Aggregates", UCAggFact22);
    coarseMapFact22->SetFactory("Nullspace", nspFact22);
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
    M22->SetFactory("UnAmalgamationInfo", amalgFact22);
    M22->SetIgnoreUserData(true);  // always use data from factories defined in factory manager

    ///////////////////////////////////////////////////////////////////////
    // define block transfer operators
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<BlockedPFactory> PFact =
        Teuchos::rcp(new BlockedPFactory());  // use row map index base from bOp
    PFact->AddFactoryManager(M11);
    PFact->AddFactoryManager(M22);

    Teuchos::RCP<GenericRFactory> RFact = Teuchos::rcp(new GenericRFactory());
    RFact->SetFactory("P", PFact);

    ///////////////////////////////////////////////////////////////////////
    // define RAPFactory
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<BlockedRAPFactory> AcFact = Teuchos::rcp(new BlockedRAPFactory());
    AcFact->SetFactory("A", MueLu::NoFactory::getRCP());  // check me!
    AcFact->SetFactory("P", PFact);
    AcFact->SetFactory("R", RFact);
    AcFact->SetRepairZeroDiagonal(true);  // repair zero diagonal entries in Ac, that are resulting
                                          // from Ptent with nullspacedim > ndofspernode

    ///////////////////////////////////////////////////////////////////////
    // transfer "SlaveDofMap" to next coarser level
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<MapTransferFactory> cmTransFact3 =
        Teuchos::rcp(new MapTransferFactory(/*"SlaveDofMap", MueLu::NoFactory::getRCP()*/));
    cmTransFact3->SetParameter("map: factory", Teuchos::ParameterEntry(std::string("NoFactory")));
    cmTransFact3->SetParameter("map: name", Teuchos::ParameterEntry(std::string("SlaveDofMap")));
    cmTransFact3->SetFactory("P", Ptent11Fact);
    AcFact->AddTransferFactory(cmTransFact3);

    ///////////////////////////////////////////////////////////////////////
    // introduce rebalancing operators
    ///////////////////////////////////////////////////////////////////////
#ifdef HAVE_MUELU_ISORROPIA

    Teuchos::RCP<RebalanceBlockAcFactory> RebalancedAcFact = Teuchos::null;
    Teuchos::RCP<RebalanceBlockInterpolationFactory> RebalancedBlockPFact = Teuchos::null;
    Teuchos::RCP<RebalanceBlockRestrictionFactory> RebalancedBlockRFact = Teuchos::null;
    if (bDoRepartition == true)
    {
      // extract subblocks from coarse unbalanced matrix

      Teuchos::RCP<SubBlockAFactory> rebA11Fact = Teuchos::rcp(new SubBlockAFactory());
      Teuchos::RCP<SubBlockAFactory> rebA22Fact = Teuchos::rcp(new SubBlockAFactory());

      rebA11Fact->SetFactory("A", AcFact);
      rebA11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
      rebA11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
      rebA22Fact->SetFactory("A", AcFact);
      rebA22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
      rebA22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

      // define rebalancing factory for coarse block matrix A(1,1)
      Teuchos::RCP<AmalgamationFactory> rebAmalgFact11 = Teuchos::rcp(new AmalgamationFactory());
      rebAmalgFact11->SetFactory("A", rebA11Fact);
      Teuchos::RCP<MueLu::IsorropiaInterface<LO, GO, NO>> isoInterface1 =
          Teuchos::rcp(new MueLu::IsorropiaInterface<LO, GO, NO>());
      isoInterface1->SetFactory("A", rebA11Fact);
      isoInterface1->SetFactory("UnAmalgamationInfo", rebAmalgFact11);
      Teuchos::RCP<MueLu::RepartitionInterface<LO, GO, NO>> repInterface1 =
          Teuchos::rcp(new MueLu::RepartitionInterface<LO, GO, NO>());
      repInterface1->SetFactory("A", rebA11Fact);
      repInterface1->SetFactory("AmalgamatedPartition", isoInterface1);

      // Repartitioning (creates "Importer" from "Partition")
      Teuchos::RCP<Factory> RepartitionFact1 = Teuchos::rcp(new RepartitionFactory());
      {
        Teuchos::ParameterList paramList;
        paramList.set("minRowsPerProcessor", optMinRowsPerProc);
        paramList.set("nonzeroImbalance", optNnzImbalance);
        paramList.set("startLevel", 1);
        RepartitionFact1->SetParameterList(paramList);
      }
      RepartitionFact1->SetFactory("A", rebA11Fact);
      RepartitionFact1->SetFactory("Partition", repInterface1);

      // define rebalancing factory for coarse block matrix A(1,1)
      Teuchos::RCP<AmalgamationFactory> rebAmalgFact22 = Teuchos::rcp(new AmalgamationFactory());
      rebAmalgFact22->SetFactory("A", rebA22Fact);


      Teuchos::RCP<MueLu::ContactSPRepartitionInterface<LO, GO, NO>> repInterface2 =
          Teuchos::rcp(new MueLu::ContactSPRepartitionInterface<LO, GO, NO>());
      repInterface2->SetFactory("A", rebA22Fact);  // use blocked matrix A as input
      repInterface2->SetFactory("AmalgamatedPartition", isoInterface1);

      // second repartition factory
      Teuchos::RCP<Factory> RepartitionFact2 = Teuchos::rcp(new RepartitionFactory());
      {
        Teuchos::ParameterList paramList;
        paramList.set("minRowsPerProcessor", 1);  // turn off repartitioning
        paramList.set("nonzeroImbalance", 1.0);
        paramList.set("startLevel", 10000);

        RepartitionFact2->SetParameterList(paramList);
      }
      RepartitionFact2->SetFactory("A", rebA22Fact);
      RepartitionFact2->SetFactory("Partition", repInterface2);
      RepartitionFact2->SetFactory("number of partitions",
          RepartitionFact1);  // use the same number of partitions as repart fact 1

      ///////////////////////////////////////////////////////////////////////
      // define rebalanced factory managers
      ///////////////////////////////////////////////////////////////////////

      Teuchos::RCP<FactoryManager> rebM11 = Teuchos::rcp(new FactoryManager());
      rebM11->SetFactory("A", AcFact);  // coarse level non-rebalanced block matrix
      rebM11->SetFactory("Importer", RepartitionFact1);
      rebM11->SetFactory("Nullspace", nspFact11);  // needed for rebalancing

      Teuchos::RCP<FactoryManager> rebM22 = Teuchos::rcp(new FactoryManager());
      rebM22->SetFactory("A", AcFact);  // coarse level non-rebalanced block matrix
      rebM22->SetFactory("Importer", RepartitionFact2);
      rebM22->SetFactory("Nullspace", nspFact22);

      // reorder transfer operators
      RebalancedBlockPFact = Teuchos::rcp(new RebalanceBlockInterpolationFactory());
      RebalancedBlockPFact->SetFactory("P", PFact);
      RebalancedBlockPFact->AddFactoryManager(rebM11);
      RebalancedBlockPFact->AddFactoryManager(rebM22);

      RebalancedBlockRFact = Teuchos::rcp(new RebalanceBlockRestrictionFactory());
      RebalancedBlockRFact->SetFactory("R", RFact);
      RebalancedBlockRFact->AddFactoryManager(rebM11);
      RebalancedBlockRFact->AddFactoryManager(rebM22);

      // rebalanced coarse level matrix
      RebalancedAcFact = Teuchos::rcp(new RebalanceBlockAcFactory());
      RebalancedAcFact->SetFactory("A", AcFact);
      RebalancedAcFact->AddFactoryManager(rebM11);
      RebalancedAcFact->AddFactoryManager(rebM22);

      // rebalance slave dof map
      Teuchos::RCP<RebalanceMapFactory> rebFact = Teuchos::rcp(new RebalanceMapFactory());
      rebFact->SetParameter("Map name", Teuchos::ParameterEntry(std::string("SlaveDofMap")));
      rebFact->SetFactory("Importer", RepartitionFact1);
      RebalancedAcFact->AddRebalanceFactory(rebFact);
    }  // end if doRepartitioning
#endif

    ///////////////////////////////////////////////////////////////////////
    // create block smoother for coarsest level
    ///////////////////////////////////////////////////////////////////////
    Teuchos::RCP<SmootherFactory> SmooFactCoarsest = GetCoarsestBlockSmootherFactory(mllist_);

    ///////////////////////////////////////////////////////////////////////
    // prepare factory managers
    ///////////////////////////////////////////////////////////////////////

    bool bIsLastLevel = false;
    std::vector<Teuchos::RCP<FactoryManager>> vecManager(maxLevels);
    for (int i = 0; i < maxLevels; i++)
    {
      Teuchos::ParameterList pp(mllist_);

      // fine/intermediate level smoother
      Teuchos::RCP<SmootherFactory> SmooFactFine = GetBlockSmootherFactory(pp, i);

      vecManager[i] = Teuchos::rcp(new FactoryManager());
      if (SmooFactFine != Teuchos::null)
        vecManager[i]->SetFactory("Smoother",
            SmooFactFine);  // Hierarchy.Setup uses TOPSmootherFactory, that only needs "Smoother"
      vecManager[i]->SetFactory("CoarseSolver", SmooFactCoarsest);
#ifdef HAVE_MUELU_ISORROPIA
      if (bDoRepartition)
      {
        // Same RAP, prolongator, and restrictor factories on all levels
        vecManager[i]->SetFactory("A", RebalancedAcFact);
        vecManager[i]->SetFactory("P", RebalancedBlockPFact);
        vecManager[i]->SetFactory("R", RebalancedBlockRFact);
      }
      else
      {
        // Same RAP, prolongator, and restrictor factories on all levels
        vecManager[i]->SetFactory("A", AcFact);
        vecManager[i]->SetFactory("P", PFact);
        vecManager[i]->SetFactory("R", RFact);
      }
#else
      // Same RAP, prolongator, and restrictor factories on all levels
      vecManager[i]->SetFactory("A", AcFact);
      vecManager[i]->SetFactory("P", PFact);
      vecManager[i]->SetFactory("R", RFact);
#endif
    }

    // use new Hierarchy::Setup routine
    if (maxLevels == 1)
    {
      bIsLastLevel =
          H->Setup(0, Teuchos::null, vecManager[0], Teuchos::null);  // 1 level "multigrid" method
    }
    else
    {
      bIsLastLevel =
          H->Setup(0, Teuchos::null, vecManager[0], vecManager[1]);  // first (finest) level
      for (int i = 1; i < maxLevels - 1; i++)
      {  // intermedium levels
        if (bIsLastLevel == true) break;
        bIsLastLevel = H->Setup(i, vecManager[i - 1], vecManager[i], vecManager[i + 1]);
      }
      if (bIsLastLevel == false)
      {  // coarsest level
        bIsLastLevel = H->Setup(
            maxLevels - 1, vecManager[maxLevels - 2], vecManager[maxLevels - 1], Teuchos::null);
      }
    }

#if 0
    { // some debug output
      // print out content of levels
      std::cout << "FINAL CONTENT of multigrid levels" << std::endl;
      Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      for(LO l = 0; l < H->GetNumLevels(); l++) {
        Teuchos::RCP<Level> coarseLevel = H->GetLevel(l);
        coarseLevel->print(*out);
      }
      std::cout << "END FINAL CONTENT of multigrid levels" << std::endl;
    } // end debug output
    Write(H);
#endif

    // set multigrid preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

    // store multigrid hierarchy
    H_ = H;
  }
  else
  {
    /* Reuse existing multigrid hierarchy
     *
     * We use the same hierarchy, but just set the new fine level matrix.
     */

    H_->setlib(Xpetra::UseEpetra);  // not very nice, but safe.
    H_->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));

    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H_));
  }

  return;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactSpPreconditioner::GetBlockSmootherFactory(
    const Teuchos::ParameterList& paramList, int level) const
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

  return CreateBlockSmootherFactory(type, smolevelsublist, false);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactSpPreconditioner::GetCoarsestBlockSmootherFactory(
    const Teuchos::ParameterList& paramList) const
{
  std::string type = paramList.get<std::string>("coarse: type");
  TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), MueLu::Exceptions::RuntimeError,
      "MueLu::ContactSpPreconditioner: no ML smoother type for coarsest level. error.");

  return CreateBlockSmootherFactory(type, paramList, true);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactSpPreconditioner::CreateBlockSmootherFactory(
    const std::string& type, const Teuchos::ParameterList& paramList, const bool bCoarse) const
{
  if (type == "SIMPLE" || type == "SIMPLEC")
  {  // SIMPLE variant
    return GetSIMPLESmootherFactory(paramList, bCoarse);
  }
  else if (type == "Braess-Sarazin")
  {  // Braess-Sarazin smoother
    return GetBraessSarazinSmootherFactory(paramList, bCoarse);
  }
  else if (type == "IBD")
  {  // indefinite blocked diagonal preconditioner
    return GetIndefBlockedDiagonalSmootherFactory(paramList, bCoarse);
  }
  else if (type == "Uzawa")
  {  // Uzawa smoother
    return GetUzawaSmootherFactory(paramList, bCoarse);
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
        "MueLu::ContactSPPreconditioner: Please set the ML_SMOOTHERMED and ML_SMOOTHERFINE "
        "parameters to SIMPLE(C) or BS in your dat file. Other smoother options are not accepted. "
        "\n Note: In fact we're using only the ML_DAMPFINE, ML_DAMPMED, ML_DAMPCOARSE as well as "
        "the ML_SMOTIMES parameters for Braess-Sarazin.");
  }

  return Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactSpPreconditioner::GetUzawaSmootherFactory(
    const Teuchos::ParameterList& paramList, const bool bCoarse) const
{
  std::string strCoarse = "coarse";
  if (bCoarse == false) strCoarse = "smoother";

  // extract type: should be Uzawa
  std::string type = paramList.get<std::string>(strCoarse + ": type");

  // extract parameters for Uzawa
  Scalar omega = paramList.get<double>(strCoarse + ": damping factor");
  int sweeps = paramList.get<int>(strCoarse + ": sweeps");

  // define indefinite blocked diagonal smoother
  Teuchos::RCP<UzawaSmoother> smootherPrototype = Teuchos::rcp(new UzawaSmoother());
  smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(sweeps));
  smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(omega));

  // define prediction smoother/solver
  Teuchos::RCP<SubBlockAFactory> A00Fact = Teuchos::rcp(new SubBlockAFactory());
  A00Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A00Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A00Fact->SetParameter("block col", Teuchos::ParameterEntry(0));

  Teuchos::RCP<SmootherPrototype> smoProtoPred = Teuchos::null;

  const Teuchos::ParameterList& PredList = paramList.sublist(strCoarse + ": Predictor list");
  std::string PredPermstrategy = "none";
  if (PredList.isSublist("Aztec Parameters") &&
      PredList.sublist("Aztec Parameters").isParameter("permutation strategy"))
    PredPermstrategy =
        PredList.sublist("Aztec Parameters").get<std::string>("permutation strategy");
  Teuchos::RCP<FactoryManager> MBpred = Teuchos::rcp(new FactoryManager());
  Teuchos::RCP<SmootherFactory> SmooPredFact = Teuchos::null;
  if (PredPermstrategy == "none")
  {
    SmooPredFact = InterpretBACIList2MueLuSmoother(PredList, A00Fact);
    MBpred->SetFactory("A", A00Fact);
  }
  else
  {
    // define permutation factory for SchurComplement equation
    Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new PermutationFactory());
    PermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
    PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);
    PermFact->SetParameter("PermutationStrategy",
        Teuchos::ParameterEntry(std::string("Local")));  // TODO Local + stridedMaps!!!
    PermFact->SetFactory("A", A00Fact);                  // use (0,0) block matrix as input for A

    SmooPredFact = InterpretBACIList2MueLuSmoother(PredList, PermFact);
    MBpred->SetFactory("A", PermFact);
  }

  MBpred->SetFactory("Smoother", SmooPredFact);  // solver for SchurComplement equation
  MBpred->SetIgnoreUserData(true);
  // smootherPrototype->SetVelocityPredictionFactoryManager(MBpred);
  smootherPrototype->AddFactoryManager(MBpred, 0);

  // create SchurComp factory (SchurComplement smoother is provided by local FactoryManager)
  Teuchos::RCP<MueLu::SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> SFact =
      Teuchos::rcp(new MueLu::SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
  SFact->SetParameter(
      "omega", Teuchos::ParameterEntry(1.0 /*omega*/));  // hard code Schur complement
  SFact->SetParameter("lumping",
      Teuchos::ParameterEntry(false /*bSimpleC*/));  // no lumping in Schur Complement approximation
  SFact->SetFactory(
      "A", MueLu::NoFactory::getRCP());  // explicitly set AFact as input factory for
                                         // SchurComplement (must be the blocked 2x2 operator)

  // define SchurComplement solver
  Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::null;
  const Teuchos::ParameterList& SchurCompList = paramList.sublist(strCoarse + ": SchurComp list");
  std::string SchurCompPermstrategy = "none";
  if (SchurCompList.isSublist("Aztec Parameters") &&
      SchurCompList.sublist("Aztec Parameters").isParameter("permutation strategy"))
    SchurCompPermstrategy =
        SchurCompList.sublist("Aztec Parameters").get<std::string>("permutation strategy");

  // setup local factory manager for SchurComplementFactory
  Teuchos::RCP<FactoryManager> MB = Teuchos::rcp(new FactoryManager());
  Teuchos::RCP<SmootherFactory> SmooSCFact = Teuchos::null;
  if (SchurCompPermstrategy == "none")
  {
    SmooSCFact = InterpretBACIList2MueLuSmoother(SchurCompList, SFact);
    MB->SetFactory("A", SFact);  // SchurCompFactory as generating factory for SchurComp equation
  }
  else
  {
    // define permutation factory for SchurComplement equation
    Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new PermutationFactory());
    PermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
    PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);
    PermFact->SetParameter("PermutationStrategy",
        Teuchos::ParameterEntry(std::string("Local")));  // TODO Local + stridedMaps!!!
    PermFact->SetFactory("A", SFact);  // use SchurComplement matrix as input for A

    SmooSCFact = InterpretBACIList2MueLuSmoother(
        SchurCompList, PermFact);   // with permutation (use permuted A)
    MB->SetFactory("A", PermFact);  // with permutation (use permuted SchurComplement matrix as A)
  }

  MB->SetFactory("Smoother", SmooSCFact);  // solver for SchurComplement equation
  MB->SetIgnoreUserData(true);
  // smootherPrototype->SetSchurCompFactoryManager(MB);  // add SC smoother information
  smootherPrototype->AddFactoryManager(MB, 1);
  smootherPrototype->SetFactory("A", MueLu::NoFactory::getRCP());

  // create smoother factory
  Teuchos::RCP<SmootherFactory> SmooFact;
  SmooFact = Teuchos::rcp(new SmootherFactory(smootherPrototype));

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if (paramList.isParameter(strCoarse + ": pre or post"))
    preorpost = paramList.get<std::string>(strCoarse + ": pre or post");

  if (preorpost == "pre")
  {
    SmooFact->SetSmootherPrototypes(smootherPrototype, Teuchos::null);
  }
  else if (preorpost == "post")
  {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smootherPrototype);
  }

  return SmooFact;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactSpPreconditioner::GetIndefBlockedDiagonalSmootherFactory(
    const Teuchos::ParameterList& paramList, const bool bCoarse) const
{
  std::string strCoarse = "coarse";
  if (bCoarse == false) strCoarse = "smoother";

  // extract type: should be IBD (indefinite blocked diagonal smoother)
  std::string type = paramList.get<std::string>(strCoarse + ": type");

  // extract parameters for IBD
  Scalar omega = paramList.get<double>(strCoarse + ": damping factor");
  int sweeps = paramList.get<int>(strCoarse + ": sweeps");

  // define indefinite blocked diagonal smoother
  Teuchos::RCP<IndefBlockedDiagonalSmoother> smootherPrototype =
      Teuchos::rcp(new IndefBlockedDiagonalSmoother());
  smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(sweeps));
  smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(omega));

  // define prediction smoother/solver
  Teuchos::RCP<SubBlockAFactory> A00Fact = Teuchos::rcp(new SubBlockAFactory());
  A00Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A00Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A00Fact->SetParameter("block col", Teuchos::ParameterEntry(0));

  Teuchos::RCP<SmootherPrototype> smoProtoPred = Teuchos::null;

  const Teuchos::ParameterList& PredList = paramList.sublist(strCoarse + ": Predictor list");
  std::string PredPermstrategy = "none";
  if (PredList.isSublist("Aztec Parameters") &&
      PredList.sublist("Aztec Parameters").isParameter("permutation strategy"))
    PredPermstrategy =
        PredList.sublist("Aztec Parameters").get<std::string>("permutation strategy");
  Teuchos::RCP<FactoryManager> MBpred = Teuchos::rcp(new FactoryManager());
  Teuchos::RCP<SmootherFactory> SmooPredFact = Teuchos::null;
  if (PredPermstrategy == "none")
  {
    SmooPredFact = InterpretBACIList2MueLuSmoother(PredList, A00Fact);
    MBpred->SetFactory("A", A00Fact);
  }
  else
  {
    // define permutation factory for SchurComplement equation
    Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new PermutationFactory());
    PermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
    PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);
    PermFact->SetParameter("PermutationStrategy",
        Teuchos::ParameterEntry(std::string("Local")));  // TODO Local + stridedMaps!!!
    PermFact->SetFactory("A", A00Fact);                  // use (0,0) block matrix as input for A

    SmooPredFact = InterpretBACIList2MueLuSmoother(PredList, PermFact);
    MBpred->SetFactory("A", PermFact);
  }

  MBpred->SetFactory("Smoother", SmooPredFact);  // solver for SchurComplement equation
  MBpred->SetIgnoreUserData(true);
  // smootherPrototype->SetVelocityPredictionFactoryManager(MBpred);
  smootherPrototype->AddFactoryManager(MBpred, 0);

  // create SchurComp factory (SchurComplement smoother is provided by local FactoryManager)
  Teuchos::RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
  SFact->SetParameter(
      "omega", Teuchos::ParameterEntry(1.0 /*omega*/));  // hard code Schur complement
  SFact->SetParameter("lumping",
      Teuchos::ParameterEntry(false /*bSimpleC*/));  // no lumping in Schur Complement approximation
  SFact->SetFactory(
      "A", MueLu::NoFactory::getRCP());  // explicitly set AFact as input factory for
                                         // SchurComplement (must be the blocked 2x2 operator)

  // define SchurComplement solver
  Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::null;
  const Teuchos::ParameterList& SchurCompList = paramList.sublist(strCoarse + ": SchurComp list");
  std::string SchurCompPermstrategy = "none";
  if (SchurCompList.isSublist("Aztec Parameters") &&
      SchurCompList.sublist("Aztec Parameters").isParameter("permutation strategy"))
    SchurCompPermstrategy =
        SchurCompList.sublist("Aztec Parameters").get<std::string>("permutation strategy");

  // setup local factory manager for SchurComplementFactory
  Teuchos::RCP<FactoryManager> MB = Teuchos::rcp(new FactoryManager());
  Teuchos::RCP<SmootherFactory> SmooSCFact = Teuchos::null;
  if (SchurCompPermstrategy == "none")
  {
    SmooSCFact = InterpretBACIList2MueLuSmoother(SchurCompList, SFact);
    MB->SetFactory("A", SFact);  // SchurCompFactory as generating factory for SchurComp equation
  }
  else
  {
    // define permutation factory for SchurComplement equation
    Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new PermutationFactory());
    PermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
    PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);
    PermFact->SetParameter("PermutationStrategy",
        Teuchos::ParameterEntry(std::string("Local")));  // TODO Local + stridedMaps!!!
    PermFact->SetFactory("A", SFact);  // use SchurComplement matrix as input for A

    SmooSCFact = InterpretBACIList2MueLuSmoother(
        SchurCompList, PermFact);   // with permutation (use permuted A)
    MB->SetFactory("A", PermFact);  // with permutation (use permuted SchurComplement matrix as A)
  }

  MB->SetFactory("Smoother", SmooSCFact);  // solver for SchurComplement equation
  MB->SetIgnoreUserData(true);
  // smootherPrototype->SetSchurCompFactoryManager(MB);  // add SC smoother information
  smootherPrototype->AddFactoryManager(MB, 1);
  smootherPrototype->SetFactory("A", MueLu::NoFactory::getRCP());

  // create smoother factory
  Teuchos::RCP<SmootherFactory> SmooFact;
  SmooFact = Teuchos::rcp(new SmootherFactory(smootherPrototype));

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if (paramList.isParameter(strCoarse + ": pre or post"))
    preorpost = paramList.get<std::string>(strCoarse + ": pre or post");

  if (preorpost == "pre")
  {
    SmooFact->SetSmootherPrototypes(smootherPrototype, Teuchos::null);
  }
  else if (preorpost == "post")
  {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smootherPrototype);
  }

  return SmooFact;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactSpPreconditioner::GetSIMPLESmootherFactory(
    const Teuchos::ParameterList& paramList, const bool bCoarse) const
{
  std::string strCoarse = "coarse";
  if (bCoarse == false) strCoarse = "smoother";

  // check whether to use SIMPLEC or SIMPLE
  bool bSimpleC = false;
  std::string type = paramList.get<std::string>(strCoarse + ": type");
  if (type == "SIMPLEC") bSimpleC = true;

  // extract parameters for SIMPLE/C
  Scalar omega = paramList.get<double>(strCoarse + ": damping factor");
  int sweeps = paramList.get<int>(strCoarse + ": sweeps");

  // define SIMPLE smoother
  Teuchos::RCP<SimpleSmoother> smootherPrototype = Teuchos::rcp(new SimpleSmoother());
  smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(sweeps));
  smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(omega));
  smootherPrototype->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(bSimpleC));

  // define prediction smoother/solver
  Teuchos::RCP<SubBlockAFactory> A00Fact = Teuchos::rcp(new SubBlockAFactory());
  A00Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A00Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A00Fact->SetParameter("block col", Teuchos::ParameterEntry(0));

  Teuchos::RCP<SmootherPrototype> smoProtoPred = Teuchos::null;

  const Teuchos::ParameterList& PredList = paramList.sublist(strCoarse + ": Predictor list");
  std::string PredPermstrategy = "none";
  if (PredList.isSublist("Aztec Parameters") &&
      PredList.sublist("Aztec Parameters").isParameter("permutation strategy"))
    PredPermstrategy =
        PredList.sublist("Aztec Parameters").get<std::string>("permutation strategy");
  Teuchos::RCP<FactoryManager> MBpred = Teuchos::rcp(new FactoryManager());
  Teuchos::RCP<SmootherFactory> SmooPredFact = Teuchos::null;
  if (PredPermstrategy == "none")
  {
    SmooPredFact = InterpretBACIList2MueLuSmoother(PredList, A00Fact);
    MBpred->SetFactory("A", A00Fact);
  }
  else
  {
    // define permutation factory for SchurComplement equation
    Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new PermutationFactory());
    PermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
    PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);
    PermFact->SetParameter("PermutationStrategy", Teuchos::ParameterEntry(std::string("Local")));
    PermFact->SetFactory("A", A00Fact);  // use (0,0) block matrix as input for A

    SmooPredFact = InterpretBACIList2MueLuSmoother(PredList, PermFact);
    MBpred->SetFactory("A", PermFact);
  }

  MBpred->SetFactory("Smoother", SmooPredFact);  // solver for SchurComplement equation
  MBpred->SetIgnoreUserData(true);
  smootherPrototype->SetVelocityPredictionFactoryManager(MBpred);

  // create SchurComp factory (SchurComplement smoother is provided by local FactoryManager)
  Teuchos::RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
  SFact->SetParameter("omega", Teuchos::ParameterEntry(omega));
  SFact->SetParameter("lumping", Teuchos::ParameterEntry(bSimpleC));
  SFact->SetFactory(
      "A", MueLu::NoFactory::getRCP());  // explicitly set AFact as input factory for
                                         // SchurComplement (must be the blocked 2x2 operator)

  // define SchurComplement solver
  Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::null;
  const Teuchos::ParameterList& SchurCompList = paramList.sublist(strCoarse + ": SchurComp list");
  std::string SchurCompPermstrategy = "none";
  if (SchurCompList.isSublist("Aztec Parameters") &&
      SchurCompList.sublist("Aztec Parameters").isParameter("permutation strategy"))
    SchurCompPermstrategy =
        SchurCompList.sublist("Aztec Parameters").get<std::string>("permutation strategy");

  // setup local factory manager for SchurComplementFactory
  Teuchos::RCP<FactoryManager> MB = Teuchos::rcp(new FactoryManager());
  Teuchos::RCP<SmootherFactory> SmooSCFact = Teuchos::null;
  if (SchurCompPermstrategy == "none")
  {
    SmooSCFact = InterpretBACIList2MueLuSmoother(SchurCompList, SFact);
    MB->SetFactory("A", SFact);  // SchurCompFactory as generating factory for SchurComp equation
  }
  else
  {
    // define permutation factory for SchurComplement equation
    Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new PermutationFactory());
    PermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
    PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);
    PermFact->SetParameter("PermutationStrategy", Teuchos::ParameterEntry(std::string("Local")));
    PermFact->SetFactory("A", SFact);  // use SchurComplement matrix as input for A

    SmooSCFact = InterpretBACIList2MueLuSmoother(
        SchurCompList, PermFact);   // with permutation (use permuted A)
    MB->SetFactory("A", PermFact);  // with permutation (use permuted SchurComplement matrix as A)
  }

  MB->SetFactory("Smoother", SmooSCFact);  // solver for SchurComplement equation
  MB->SetIgnoreUserData(true);
  smootherPrototype->SetSchurCompFactoryManager(MB);  // add SC smoother information
  smootherPrototype->SetFactory("A", MueLu::NoFactory::getRCP());

  // create smoother factory
  Teuchos::RCP<SmootherFactory> SmooFact;
  SmooFact = Teuchos::rcp(new SmootherFactory(smootherPrototype));

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if (paramList.isParameter(strCoarse + ": pre or post"))
    preorpost = paramList.get<std::string>(strCoarse + ": pre or post");

  if (preorpost == "pre")
  {
    SmooFact->SetSmootherPrototypes(smootherPrototype, Teuchos::null);
  }
  else if (preorpost == "post")
  {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smootherPrototype);
  }

  return SmooFact;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactSpPreconditioner::GetBraessSarazinSmootherFactory(
    const Teuchos::ParameterList& paramList, const bool bCoarse) const
{
  std::string strCoarse = "coarse";
  if (bCoarse == false) strCoarse = "smoother";

  std::string type = paramList.get<std::string>(strCoarse + ": type");
  TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), MueLu::Exceptions::RuntimeError,
      "MueLu::Interpreter: no ML smoother type for level. error.");

  if (type != "Braess-Sarazin")
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
        "MueLu::ContactSPPreconditioner: Please set the ML_SMOOTHERCOARSE, ML_SMOOTHERMED and "
        "ML_SMOOTHERFINE parameters to BS in your dat file. Other smoother options are not "
        "accepted. \n Note: In fact we're using only the ML_DAMPFINE, ML_DAMPMED, ML_DAMPCOARSE as "
        "well as the ML_SMOTIMES parameters for Braess-Sarazin.");
  }

  // extract parameters for Braess-Sarazin smoother
  Scalar omega = paramList.get<double>(
      strCoarse + ": damping factor");  // Braess-Sarazin damping/scaling factor
  int sweeps = paramList.get<int>(strCoarse + ": sweeps");

  // create SchurComp factory (SchurComplement smoother is provided by local FactoryManager)
  Teuchos::RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
  SFact->SetParameter("omega", Teuchos::ParameterEntry(omega));
  SFact->SetParameter("lumping", Teuchos::ParameterEntry(false));
  SFact->SetFactory("A", MueLu::NoFactory::getRCP()); /* XXX check me */
  Teuchos::RCP<BraessSarazinSmoother> smootherPrototype =
      Teuchos::rcp(new BraessSarazinSmoother());  // append SC smoother information
  smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(sweeps));
  smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(omega));

  // define SchurComplement solver
  const Teuchos::ParameterList& SchurCompList = paramList.sublist(strCoarse + ": SchurComp list");
  std::string permstrategy = "none";
  if (SchurCompList.isSublist("Aztec Parameters") &&
      SchurCompList.sublist("Aztec Parameters").isParameter("permutation strategy"))
    permstrategy =
        SchurCompList.sublist("Aztec Parameters").get<std::string>("permutation strategy");


  // define SchurComplement smoother and
  // setup local factory manager for SchurComplementFactory
  Teuchos::RCP<SmootherFactory> SmooSCFact = Teuchos::null;
  Teuchos::RCP<FactoryManager> MB = Teuchos::rcp(new FactoryManager());
  if (permstrategy == "none")
  {
    SmooSCFact = InterpretBACIList2MueLuSmoother(SchurCompList, SFact);
    MB->SetFactory("A", SFact);  // SchurCompFactory as generating factory for SchurComp equation
  }
  else
  {
    // define permutation factory for SchurComplement equation
    Teuchos::RCP<PermutationFactory> PermFact = Teuchos::rcp(new PermutationFactory());
    PermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
    PermFact->SetFactory("PermutationRowMapFactory", Teuchos::null);
    PermFact->SetParameter("PermutationStrategy",
        Teuchos::ParameterEntry(std::string("Local")));  // TODO Local + stridedMaps!!!
    PermFact->SetFactory("A", SFact);  // use SchurComplement matrix as input for A

    SmooSCFact = InterpretBACIList2MueLuSmoother(
        SchurCompList, PermFact);   // with permutation (use permuted A)
    MB->SetFactory("A", PermFact);  // with permutation (use permuted SchurComplement matrix as A)
  }

  // share common code
  MB->SetFactory("Smoother", SmooSCFact);  // solver for SchurComplement equation
  MB->SetIgnoreUserData(true);
  smootherPrototype->AddFactoryManager(MB, 0);  // add SC smoother information
  smootherPrototype->SetFactory("A", MueLu::NoFactory::getRCP());

  // create smoother factory
  Teuchos::RCP<SmootherFactory> SmooFact;
  SmooFact = Teuchos::rcp(new SmootherFactory(smootherPrototype));

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if (paramList.isParameter(strCoarse + ": pre or post"))
    preorpost = paramList.get<std::string>(strCoarse + ": pre or post");

  if (preorpost == "pre")
  {
    SmooFact->SetSmootherPrototypes(smootherPrototype, Teuchos::null);
  }
  else if (preorpost == "post")
  {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smootherPrototype);
  }

  return SmooFact;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
LINALG::SOLVER::MueLuContactSpPreconditioner::InterpretBACIList2MueLuSmoother(
    const Teuchos::ParameterList& paramList, Teuchos::RCP<FactoryBase> AFact) const
{
  Teuchos::RCP<SmootherPrototype> smoProtoSC = Teuchos::null;
  std::string type = paramList.get<std::string>("solver");
  if (type == "umfpack")
  {
    smoProtoSC = Teuchos::rcp(new DirectSolver("Klu" /*"Umfpack"*/, Teuchos::ParameterList()));
  }
  else if (type == "aztec")
  {
    std::string prectype =
        paramList.sublist("Aztec Parameters").get<std::string>("Preconditioner Type");
    std::string permstrategy = "none";
    if (paramList.isSublist("Aztec Parameters") &&
        paramList.sublist("Aztec Parameters").isParameter("permutation strategy"))
      permstrategy = paramList.sublist("Aztec Parameters").get<std::string>("permutation strategy");
    if (permstrategy == "none")
    {
      // no permutation strategy: build plain smoothers
      if (prectype == "ILU")
      {
        Teuchos::ParameterList SCList;
        SCList.set<int>("fact: level-of-fill", 0);
        smoProtoSC =
            MueLu::GetIfpackSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>("ILU", SCList, 0);
        smoProtoSC->SetFactory("A", AFact);
      }
      else if (prectype == "point relaxation")
      {
        const Teuchos::ParameterList& ifpackParameters = paramList.sublist("IFPACK Parameters");
        smoProtoSC = Teuchos::rcp(new TrilinosSmoother("RELAXATION", ifpackParameters, 0));
        smoProtoSC->SetFactory("A", AFact);
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
            "MueLu::ContactSPPreconditioner: unknown type for SchurComplement smoother (IFPACK "
            "based)");
    }
    else
    {
      // use PermutingSmoother
      if (prectype == "ILU")
      {
        const Teuchos::ParameterList& ifpackParameters = paramList.sublist("IFPACK Parameters");
        smoProtoSC = Teuchos::rcp(
            new PermutingSmoother("", Teuchos::null, "ILU", ifpackParameters, 0, AFact));
      }
      else if (prectype == "point relaxation")
      {
        const Teuchos::ParameterList& ifpackParameters = paramList.sublist("IFPACK Parameters");
        smoProtoSC = Teuchos::rcp(
            new PermutingSmoother("", Teuchos::null, "RELAXATION", ifpackParameters, 0, AFact));
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
            "MueLu::ContactSPPreconditioner: unknown type for SchurComplement smoother (IFPACK "
            "based)");
    }
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
        "MueLu::ContactSPPreconditioner: unknown type for SchurComplement solver/smoother");
  }
  smoProtoSC->SetFactory("A", AFact);
  Teuchos::RCP<SmootherFactory> SmooSCFact = Teuchos::rcp(new SmootherFactory(smoProtoSC));
  return SmooSCFact;
}

void LINALG::SOLVER::MueLuContactSpPreconditioner::Write(const Teuchos::RCP<Hierarchy>& H)
{
  LocalOrdinal startLevel = 0;
  LocalOrdinal endLevel = H->GetNumLevels() - 1;

  for (LO i = startLevel; i < endLevel + 1; ++i)
  {
    std::ostringstream buf;
    buf << i;
    std::string fileName = "A_" + buf.str() + ".m";

    Teuchos::RCP<Matrix> A = H->GetLevel(i)->Get<Teuchos::RCP<Matrix>>("A");
    Teuchos::RCP<BlockedCrsMatrix> Ab = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A);
    Teuchos::RCP<CrsMatrix> Am = Ab->Merge();
    Teuchos::RCP<Matrix> Aw = Teuchos::rcp(new CrsMatrixWrap(Am));

    Utils::Write(fileName, *Aw);

    if (i > 0)
    {
      fileName = "P_" + buf.str() + ".m";
      // Utils::Write( fileName,*(Levels_[i]-> template Get< Teuchos::RCP< Matrix> >("P")) );

      Teuchos::RCP<Matrix> P = H->GetLevel(i)->Get<Teuchos::RCP<Matrix>>("P");
      Teuchos::RCP<BlockedCrsMatrix> Pb = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(P);
      Teuchos::RCP<CrsMatrix> Pm = Pb->Merge();
      Teuchos::RCP<Matrix> Pw = Teuchos::rcp(new CrsMatrixWrap(Pm));
      Utils::Write(fileName, *Pw);

      // if (!implicitTranspose_) {

      fileName = "R_" + buf.str() + ".m";
      // Utils::Write( fileName,*(Levels_[i]-> template Get< Teuchos::RCP< Matrix> >("R")) );
      Teuchos::RCP<Matrix> R = H->GetLevel(i)->Get<Teuchos::RCP<Matrix>>("R");
      Teuchos::RCP<BlockedCrsMatrix> Rb = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(R);
      Teuchos::RCP<CrsMatrix> Rm = Rb->Merge();
      Teuchos::RCP<Matrix> Rw = Teuchos::rcp(new CrsMatrixWrap(Rm));
      Utils::Write(fileName, *Rw);

      //}
    }
  }

}  // Write()


#endif  // HAVE_MueLu
