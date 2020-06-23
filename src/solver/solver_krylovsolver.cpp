/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of Baci's interface to Krylov solvers

\level 0

\maintainer Martin Kronbichler
*/
/*---------------------------------------------------------------------*/

#include <MueLu_ConfigDefs.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <MueLu.hpp>
#include <MueLu_FactoryBase.hpp>
#include <MueLu_PermutationFactory.hpp>
#include <MueLu_SmootherPrototype.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_DirectSolver.hpp>
#ifdef TRILINOS_Q1_2015
#include <MueLu_HierarchyHelpers.hpp>
#else
#include <MueLu_HierarchyUtils.hpp>
#endif
#include <MueLu_VerboseObject.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixWrap;
typedef Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Vector;
typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactory;
typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
typedef Xpetra::EpetraCrsMatrix EpetraCrsMatrix;

typedef MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> FactoryManager;
typedef MueLu::FactoryBase FactoryBase;

typedef Scalar SC;
typedef LocalOrdinal LO;
typedef GlobalOrdinal GO;
typedef Node NO;

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <az_aztec_defs.h>  // for AZ_none (provokes compiler warning due to redeclaration of HAVE_SYS_TIME_H in mpi.h and AztecOO_config.h -> AztecOO problem)

#include "../drt_lib/drt_dserror.H"
#include "solver_krylovsolver.H"

#include "solver_pointpreconditioner.H"
#include "solver_blockpreconditioners.H"
#include "solver_krylovprojectionpreconditioner.H"
#include "solver_ifpackpreconditioner.H"
#include "solver_mlpreconditioner.H"
#include "solver_muelupreconditioner.H"
#include "solver_amgnxn_preconditioner.H"
#ifdef TRILINOS_Q1_2015
#include "solver_muelucontactpreconditioner.H"
#include "solver_muelucontactpreconditioner2.H"
#include "solver_muelucontactsppreconditioner.H"
#include "solver_muelucontactpenaltypreconditioner.H"
#endif  // TRILINOS_Q1_2015
#ifdef HAVE_TEKO
#include "solver_tekopreconditioner.H"
#endif  // HAVE_TEKO

#include <Teuchos_TimeMonitor.hpp>

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::KrylovSolver::KrylovSolver(
    const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile)
    : comm_(comm),
      params_(params),
      outfile_(outfile),
      ncall_(0),
      activeDofMap_(Teuchos::null),
      bAllowPermutation_(false),
      bPermuteLinearSystem_(false),
      permutationStrategy_("none"),
      diagDominanceRatio_(1.0),
      PermFact_(Teuchos::null)
{
  data_ = Teuchos::rcp(new MueLu::Level());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::KrylovSolver::~KrylovSolver()
{
  preconditioner_ = Teuchos::null;
  A_ = Teuchos::null;
  x_ = Teuchos::null;
  b_ = Teuchos::null;
  activeDofMap_ = Teuchos::null;
  data_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::KrylovSolver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)
{
  return preconditioner_->ApplyInverse(X, Y);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
bool LINALG::SOLVER::KrylovSolver::AllowReusePreconditioner(const int reuse, const bool reset)
{
  bool bAllowReuse = true;  // default: allow reuse of preconditioner

  // first, check some parameters with information that has to be updated
  Teuchos::ParameterList* linSysParams = NULL;
  if (Params().isSublist("Aztec Parameters"))
  {
    Teuchos::ParameterList& reflinSysParams = Params().sublist("Aztec Parameters");
    linSysParams = &reflinSysParams;
  }
  else if (Params().isSublist("Belos Parameters"))
  {
    Teuchos::ParameterList& reflinSysParams = Params().sublist("Belos Parameters");
    linSysParams = &reflinSysParams;
  }

  CheckReuseStatusOfActiveSet(bAllowReuse, linSysParams);

  // true, if preconditioner must not reused but is to re-created!
  const bool create = reset or not Ncall() or not reuse or (Ncall() % reuse) == 0;
  if (create) bAllowReuse = false;  // we have to create a new preconditioner

  // here, each processor has its own local decision made
  // bAllowReuse = true -> preconditioner can be reused
  // bAllowReuse = false -> preconditioner has to be recomputed
  // If one or more processors decide that the preconditioner has to be recomputed
  // all of the processors have to recompute it

  // synchronize results of all processors
  // all processors have to do the same (either recompute preconditioner or allow reusing it)
  int nProc = comm_.NumProc();
  int lAllowReuse = bAllowReuse == true ? 1 : 0;
  int gAllowReuse = 0;
  comm_.SumAll(&lAllowReuse, &gAllowReuse, 1);

  if (gAllowReuse == nProc)
    bAllowReuse = true;
  else
    bAllowReuse = false;

  return bAllowReuse;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovSolver::CheckReuseStatusOfActiveSet(
    bool& bAllowReuse, const Teuchos::ParameterList* linSysParams)
{
  if (linSysParams != NULL)
  {
    if (linSysParams->isSublist("Linear System properties"))
    {
      const Teuchos::ParameterList& linSystemProps =
          linSysParams->sublist("Linear System properties");

      if (linSystemProps.isParameter("contact activeDofMap"))
      {
        Teuchos::RCP<Epetra_Map> epActiveDofMap = Teuchos::null;
        epActiveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact activeDofMap");

        // Do we have history information available?
        if (activeDofMap_.is_null())
        {
          /* No history available.
           * This is the first application of the preconditioner. We cannot reuse it.
           */
          bAllowReuse = false;
        }
        else
        {
          /* History is available. We actually have to check for a change in the active set
           * by comparing the current map of active DOFs with the stored map of active DOFs
           * from the previous application of the preconditioner.
           */
          if (not epActiveDofMap->PointSameAs(*activeDofMap_))
          {
            // Map of active nodes has changed -> force preconditioner to be rebuilt
            bAllowReuse = false;
          }
        }

        // Store current map of active slave DOFs for comparison in next preconditioner application
        activeDofMap_ = epActiveDofMap;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovSolver::CreatePreconditioner(Teuchos::ParameterList& azlist,
    const bool isCrsMatrix, Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Solver:  1.1)   CreatePreconditioner");

  preconditioner_ = Teuchos::null;

  if (isCrsMatrix)
  {
    // get type of preconditioner and build either Ifpack or ML
    // if we have an ifpack parameter list, we do ifpack
    // if we have an ml parameter list we do ml
    // if we have a downwinding flag we downwind the linear problem
    if (Params().isSublist("IFPACK Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new LINALG::SOLVER::IFPACKPreconditioner(
          outfile_, Params().sublist("IFPACK Parameters"), azlist));
    }
    else if (Params().isSublist("ML Parameters"))
    {
      preconditioner_ = Teuchos::rcp(
          new LINALG::SOLVER::MLPreconditioner(outfile_, Params().sublist("ML Parameters")));
    }
    else if (Params().isSublist("MueLu Parameters"))
    {
      preconditioner_ = Teuchos::rcp(
          new LINALG::SOLVER::MueLuPreconditioner(outfile_, Params().sublist("MueLu Parameters")));
    }
    else if (Params().isSublist("MueLu (Contact) Parameters"))
    {
#ifdef TRILINOS_Q1_2015
      preconditioner_ = Teuchos::rcp(new LINALG::SOLVER::MueLuContactPreconditioner(
          outfile_, Params().sublist("MueLu (Contact) Parameters")));
#else
      dserror("MueLu (Contact) preconditioner only available with Trilinos Q1_2015.");
#endif
    }
    else if (Params().isSublist("MueLu (Contact2) Parameters"))
    {
#ifdef TRILINOS_Q1_2015
      preconditioner_ = Teuchos::rcp(new LINALG::SOLVER::MueLuContactPreconditioner2(
          outfile_, Params().sublist("MueLu (Contact2) Parameters")));
#else
      dserror("MueLu (Contact2) preconditioner only available with Trilinos Q1_2015.");
#endif
    }
    else if (Params().isSublist("MueLu (PenaltyContact) Parameters"))
    {
#ifdef TRILINOS_Q1_2015
      preconditioner_ = Teuchos::rcp(new LINALG::SOLVER::MueLuContactPenaltyPreconditioner(
          outfile_, Params().sublist("MueLu (PenaltyContact) Parameters")));
#else
      dserror("MueLu (PenaltyContact) preconditioner only available with Trilinos Q1_2015.");
#endif
    }
    else if (azlist.get<int>("AZ_precond") == AZ_none)  // FIXME Attention: this is dangerous.
    {
      preconditioner_ = Teuchos::rcp(new LINALG::SOLVER::NonePreconditioner(outfile_, Params()));
    }
    else
    {
      dserror("unknown preconditioner");
    }

    // decide whether we do what kind of scaling
    std::string scaling = azlist.get("scaling", "none");
    if (scaling == "none")
    {
    }
    else if (scaling == "infnorm")
    {
      preconditioner_ = Teuchos::rcp(new LINALG::SOLVER::InfNormPreconditioner(preconditioner_));
    }
    else if (scaling == "symmetric")
    {
      preconditioner_ = Teuchos::rcp(new LINALG::SOLVER::SymDiagPreconditioner(preconditioner_));
    }
    else
      dserror("Unknown type of scaling found in parameter list");

    if (azlist.get<bool>("downwinding", false))
    {
      preconditioner_ =
          Teuchos::rcp(new LINALG::SOLVER::DWindPreconditioner(outfile_, preconditioner_, azlist));
    }

    if (projector != Teuchos::null)
    {
      preconditioner_ = Teuchos::rcp(
          new LINALG::SOLVER::KrylovProjectionPreconditioner(outfile_, preconditioner_, projector));
    }
  }
  else
  {
    // assume block matrix

    if (Params().isSublist(
            "SIMPLER"))  // old BACI::(Cheap)SIMPLER preconditioner TODO: remove/replace me
    {
      dserror("SIMPLER sublist not supported any more.");
      preconditioner_ =
          Teuchos::rcp(new SimplePreconditioner(outfile_, Params()));  // Michael's SIMPLE for Fluid
    }
    else if (Params().isSublist("CheapSIMPLE Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new SimplePreconditioner(outfile_, Params()));
    }
    else if (Params().isSublist("BGS Parameters"))
    {
      preconditioner_ = Teuchos::rcp(
          new BGSPreconditioner(outfile_, Params(), Params().sublist("BGS Parameters")));
    }
    else if (Params().isSublist("Teko Parameters"))
    {
#ifdef HAVE_TEKO
      preconditioner_ = Teuchos::rcp(new TekoPreconditioner(outfile_, Params()));
#else
      dserror("You need the HAVE_TEKO define flag set. Works only for Trilinos Q1/2012 or newer.");
#endif
    }
    else if (Params().isSublist("MueLu Parameters"))
    {
      preconditioner_ = Teuchos::rcp(
          new MueLuBlockPreconditioner(outfile_, Params().sublist("MueLu Parameters")));
    }
    else if (Params().isSublist("MueLu (Contact) Parameters"))
    {
#ifdef TRILINOS_Q1_2015
      preconditioner_ = Teuchos::rcp(new LINALG::SOLVER::MueLuContactSpPreconditioner(
          outfile_, Params().sublist("MueLu (Contact) Parameters")));
#else
      dserror("MueLu (Contact) preconditioner only available with Trilinos Q1_2015.");
#endif
    }
    else if (Params().isSublist("AMGnxn Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new LINALG::SOLVER::AMGnxn_Preconditioner(outfile_, Params()));
    }
    else
    {
      dserror("unknown preconditioner for block matrix solver");
    }
  }

#if 0
  preconditioner_->Print( std::cout );
  std::cout << "\n";
#endif
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovSolver::BuildPermutationOperator(
    const Teuchos::RCP<Epetra_CrsMatrix>& A, const Teuchos::RCP<Epetra_Map>& epSlaveDofMap)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A));
  Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>> xCrsOp =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>(xCrsA));
  Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>> xOp =
      Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LO, GO, Node>>(xCrsOp);
  xOp->SetFixedBlockSize(Params()
                             .sublist("NodalBlockInformation")
                             .get<int>("number of momentum dofs"));  // set nBlockSize

  data_->setDefaultVerbLevel(Teuchos::VERB_NONE);
  data_->setlib(Xpetra::UseEpetra);
  data_->Set("A", xOp);


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // check, if "SlaveDofMap" information is available in parameter lists
  if (epSlaveDofMap != Teuchos::null)
  {
#ifdef TRILINOS_Q1_2015
    Teuchos::RCP<Xpetra::EpetraMap> xSlaveDofMap =
        Teuchos::rcp(new Xpetra::EpetraMap(epSlaveDofMap));
#else
    Teuchos::RCP<Xpetra::EpetraMapT<GO, NO>> xSlaveDofMap =
        Teuchos::rcp(new Xpetra::EpetraMapT<GO, NO>(epSlaveDofMap));
#endif
    data_->Set("SlaveDofMap", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, Node>>(
                                  xSlaveDofMap));  // set map with active dofs

    // define permutation factory for permuting the full matrix A
    PermFact_ =
        Teuchos::rcp(new MueLu::PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
    PermFact_->SetParameter(
        "PermutationRowMapName", Teuchos::ParameterEntry(std::string("SlaveDofMap")));
    PermFact_->SetFactory("PermutationRowMapFactory", MueLu::NoFactory::getRCP());
    PermFact_->SetParameter("PermutationStrategy", Teuchos::ParameterEntry(permutationStrategy_));
  }
  else
  {
    // permute full matrix
    PermFact_ =
        Teuchos::rcp(new MueLu::PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
    PermFact_->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
    PermFact_->SetFactory("PermutationRowMapFactory", Teuchos::null);
    PermFact_->SetParameter("PermutationStrategy", Teuchos::ParameterEntry(permutationStrategy_));
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // setup main factory manager
  Teuchos::RCP<FactoryManager> M = Teuchos::rcp(new FactoryManager());
  M->SetFactory("permQT", PermFact_);
  M->SetFactory("A", MueLu::NoFactory::getRCP());  // this is the input matrix
  MueLu::SetFactoryManager SFMFinest(data_, M);    // set factory manager for data container

  // prepare building process for permutation operators
  data_->Request("A", PermFact_.get());
  data_->Request("permA", PermFact_.get());
  data_->Request("permP", PermFact_.get());
  data_->Request("permQT", PermFact_.get());
  data_->Request("permScaling", PermFact_.get());
  data_->Request("#RowPermutations", PermFact_.get());
  data_->Request("#ColPermutations", PermFact_.get());
  data_->Request("#WideRangeRowPermutations", PermFact_.get());
  data_->Request("#WideRangeColPermutations", PermFact_.get());

  // build permutation operators
  PermFact_->Build(*data_);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovSolver::PermuteLinearSystem(
    const Teuchos::RCP<Epetra_CrsMatrix>& A, const Teuchos::RCP<Epetra_MultiVector>& b)
{
  if (!data_->IsAvailable("A", PermFact_.get()) || !data_->IsAvailable("permP", PermFact_.get()) ||
      !data_->IsAvailable("permScaling", PermFact_.get()))
    dserror("PermutedAztecSolver: call BuildPermutationOperator before PermuteLinearSystem");

  // extract permutation operators from data_
  Teuchos::RCP<Epetra_CrsMatrix> xEpPermCrsMat = GetOperatorNonConst("A", PermFact_);
  Teuchos::RCP<const Epetra_CrsMatrix> epPermPMatrix =
      GetOperator("permP", PermFact_);  // row permutation matrix
  Teuchos::RCP<const Epetra_CrsMatrix> epPermScalingMatrix =
      GetOperator("permScaling", PermFact_);  // leftScaling matrix

  // P_trafo*b
  Teuchos::RCP<Epetra_MultiVector> btemp1 = Teuchos::rcp(new Epetra_MultiVector(*b));
  epPermPMatrix->Multiply(false, *b, *btemp1);
  // P_scaling * P_trafo * b
  epPermScalingMatrix->Multiply(false, *btemp1, *b);

  // set
  // b_ = permP * b;
  // A_ = permQ^T * A * permP
  b_ = b;              // note b is permuted
  A_ = xEpPermCrsMat;  // xEpPermCrsMat->getEpetra_CrsMatrixNonConst();
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovSolver::PermuteNullSpace(const Teuchos::RCP<Epetra_CrsMatrix>& A)
{
  // note: we usually do not permute the null space!!!

  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xOp = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xOp->SetFixedBlockSize(Params()
                             .sublist("NodalBlockInformation")
                             .get<int>("number of momentum dofs"));  // set nBlockSize

  // detect MueLu/ML Paramter list
  std::string MultiGridParameterListName = "";
  if (Params().isSublist("ML Parameters"))
    MultiGridParameterListName = "ML Parameters";
  else if (Params().isSublist("MueLu Parameters"))
    MultiGridParameterListName = "MueLu Parameters";
  else if (Params().isSublist("MueLu (Contact) Parameters"))
    MultiGridParameterListName = "MueLu (Contact) Parameters";
  else if (Params().isSublist("MueLu (Contact2) Parameters"))
    MultiGridParameterListName = "MueLu (Contact2) Parameters";
  else if (Params().isSublist("MueLu (Contact3) Parameters"))
    MultiGridParameterListName = "MueLu (Contact3) Parameters";
  else if (Params().isSublist("MueLu (PenaltyContact) Parameters"))
    MultiGridParameterListName = "MueLu (PenaltyContact) Parameters";

  // retransform nullspace vectors
  if (MultiGridParameterListName == "") return;  // no nullspace to permute

  Teuchos::RCP<Matrix> xPermQtMatrix = data_->Get<Teuchos::RCP<Matrix>>("permQT", PermFact_.get());
  int numdf = Params().sublist(MultiGridParameterListName).get<int>("PDE equations", -1);
  int dimns = Params().sublist(MultiGridParameterListName).get<int>("null space: dimension", -1);
  if (dimns == -1 || numdf == -1)
    dserror(
        "PermutedAztecSolver: Error in MueLu/ML parameters: PDE equations or null space dimension "
        "wrong.");
  Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> rowMap = xOp->getRowMap();

  Teuchos::RCP<MultiVector> nspVector =
      Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
          rowMap, dimns, true);
  Teuchos::RCP<std::vector<double>> nsdata =
      Params()
          .sublist(MultiGridParameterListName)
          .get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);

  for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
  {
    Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
    const size_t myLength = nspVector->getLocalLength();
    for (size_t j = 0; j < myLength; j++)
    {
      nspVectori[j] = (*nsdata)[i * myLength + j];
    }
  }

  // calculate transformed nullspace multivector
  Teuchos::RCP<MultiVector> permutedNspVector =
      Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
          xOp->getDomainMap(), dimns, true);
  // calculate Q * b_f
  xPermQtMatrix->apply(*nspVector, *permutedNspVector, Teuchos::TRANS);

  // write data back
  for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
  {
    Teuchos::ArrayRCP<Scalar> permutedNspVectori = permutedNspVector->getDataNonConst(i);
    const size_t myLength = permutedNspVector->getLocalLength();
    for (size_t j = 0; j < myLength; j++)
    {
      (*nsdata)[i * myLength + j] = permutedNspVectori[j];
    }
  }
#if 0
    // experiment
    Teuchos::RCP<MultiVector> nspVector2 = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
    Teuchos::RCP<std::vector<double> > nsdata2 = Params().sublist(MultiGridParameterListName).get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

    for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
        Teuchos::ArrayRCP<Scalar> nspVector2i = nspVector2->getDataNonConst(i);
        const size_t myLength = nspVector2->getLocalLength();
        for(size_t j=0; j<myLength; j++) {
                nspVector2i[j] = (*nsdata2)[i*myLength+j];
        }
    }

    Teuchos::RCP<MultiVector> test = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xOp->getRowMap(),dimns,true);
    xPermQtMatrix->apply(*nspVector2, *test, Teuchos::NO_TRANS);
    test->update(-1.0, *nspVector, 1.0);
    std::cout << *test << std::endl;
#endif
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Map> LINALG::SOLVER::KrylovSolver::FindNonDiagonalDominantRows(
    const Teuchos::RCP<Matrix>& xA, double diagDominanceRatio)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

  const Teuchos::RCP<const Teuchos::Comm<int>> comm = xA->getRowMap()->getComm();

  Teuchos::RCP<Vector> diagAVec = VectorFactory::Build(xA->getRowMap(), true);
  xA->getLocalDiagCopy(*diagAVec);
  Teuchos::ArrayRCP<const Scalar> diagAVecData = diagAVec->getData(0);
  std::vector<GlobalOrdinal>
      NonDiagonalDominantGIDs;  // vector with nondiagonaldominant row GIDs on cur proc

  // loop over all local rows in matrix A and keep diagonal entries if corresponding
  // matrix rows are not contained in permRowMap
  for (size_t row = 0; row < xA->getRowMap()->getNodeNumElements(); row++)
  {
    GlobalOrdinal grow = xA->getRowMap()->getGlobalElement(row);

    // extract local row information from matrix
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    xA->getLocalRowView(row, indices, vals);

    // find column entry with max absolute value
    Scalar maxVal = 0.0;
    for (size_t j = 0; j < Teuchos::as<size_t>(indices.size()); j++)
    {
      if (std::abs(vals[j]) > maxVal)
      {
        maxVal = std::abs(vals[j]);
      }
    }

    // check the ratio nof the diagonal entry and the entry with
    // maximal absolute value in the current row
    // if the row is diagonal dominant the ratio is 1.0
    // if the row is not diagonal dominant, the ratio is < 1.0
    // if the row has a zero on the diagonal the ratio is zero
    // if the row is a zero row (-> singular matrix) we divide zero by zero
    if (std::abs(diagAVecData[row]) / maxVal < diagDominanceRatio)
    {  // optimal would be 1.0 -> row is diagonal dominant
      NonDiagonalDominantGIDs.push_back(grow);
    }
  }

  const Teuchos::ArrayView<const LocalOrdinal> NonDiagonalDominantGIDs_view(
      &NonDiagonalDominantGIDs[0], NonDiagonalDominantGIDs.size());

  Teuchos::RCP<Map> NonDiagonalDominantGIDsMap = MapFactory::Build(xA->getRowMap()->lib(),
      Teuchos::OrdinalTraits<int>::invalid(), NonDiagonalDominantGIDs_view, 0, comm);

  return NonDiagonalDominantGIDsMap;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Map> LINALG::SOLVER::KrylovSolver::FindNonDiagonalDominantRows(
    const Teuchos::RCP<Epetra_CrsMatrix>& A, double diagDominanceRatio)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xA = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xA->SetFixedBlockSize(Params()
                            .sublist("NodalBlockInformation")
                            .get<int>("number of momentum dofs"));  // set nBlockSize

  return FindNonDiagonalDominantRows(xA, diagDominanceRatio);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Map> LINALG::SOLVER::KrylovSolver::FindZeroDiagonalEntries(
    const Teuchos::RCP<Matrix>& xA, double tolerance)
{
  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

  // extract permuted but unscaled matrix
  // Teuchos::RCP<Matrix> xPermA = data_->Get<Teuchos::RCP<Matrix> >("permA", PermFact_.get());

  const Teuchos::RCP<const Teuchos::Comm<int>> comm = xA->getRowMap()->getComm();

  Teuchos::RCP<Vector> diagAVec = VectorFactory::Build(xA->getRowMap(), true);
  xA->getLocalDiagCopy(*diagAVec);
  Teuchos::ArrayRCP<const Scalar> diagAVecData = diagAVec->getData(0);
  LocalOrdinal lNumZeros = 0;
  std::vector<GlobalOrdinal> zeroGids;
  for (size_t i = 0; i < diagAVec->getMap()->getNodeNumElements(); ++i)
  {
    if (std::abs(diagAVecData[i]) < tolerance)
    {  // pick out all rows with very small diagonal entries
      lNumZeros++;
      zeroGids.push_back(diagAVec->getMap()->getGlobalElement(i));
    }
  }

  const Teuchos::ArrayView<const LocalOrdinal> zeroGids_view(&zeroGids[0], zeroGids.size());

  Teuchos::RCP<Map> zeroDiagonalMap = MapFactory::Build(
      xA->getRowMap()->lib(), Teuchos::OrdinalTraits<int>::invalid(), zeroGids_view, 0, comm);

  return zeroDiagonalMap;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Map> LINALG::SOLVER::KrylovSolver::FindZeroDiagonalEntries(
    const Teuchos::RCP<Epetra_CrsMatrix>& A, double tolerance)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xA = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xA->SetFixedBlockSize(Params()
                            .sublist("NodalBlockInformation")
                            .get<int>("number of momentum dofs"));  // set nBlockSize

  return FindZeroDiagonalEntries(xA, tolerance);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::KrylovSolver::CountZerosOnDiagonalEpetra(
    const Teuchos::RCP<Epetra_CrsMatrix>& A)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xOp = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xOp->SetFixedBlockSize(Params()
                             .sublist("NodalBlockInformation")
                             .get<int>("number of momentum dofs"));  // set nBlockSize

  return CountZerosOnDiagonal(xOp);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::KrylovSolver::CountZerosOnDiagonal(const Teuchos::RCP<const Matrix>& xOp)
{
  Teuchos::RCP<Vector> diagAVec = VectorFactory::Build(xOp->getRowMap(), true);
  xOp->getLocalDiagCopy(*diagAVec);
  Teuchos::ArrayRCP<const Scalar> diagAVecData = diagAVec->getData(0);
  LocalOrdinal lNumZeros = 0;
  GlobalOrdinal gNumZeros = 0;
  for (size_t i = 0; i < diagAVec->getMap()->getNodeNumElements(); ++i)
  {
    if (diagAVecData[i] == 0.0)
    {
      lNumZeros++;
    }
  }

  // sum up all entries in multipleColRequests over all processors
  MueLu_sumAll(diagAVec->getMap()->getComm(), (LocalOrdinal)lNumZeros, gNumZeros);

  return Teuchos::as<int>(gNumZeros);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovSolver::ReTransformSolution()
{
  Teuchos::RCP<const Epetra_CrsMatrix> epPermQTMatrix = GetOperator("permQT", PermFact_);
  Teuchos::RCP<Epetra_MultiVector> xtemp = Teuchos::rcp(new Epetra_MultiVector(*x_));
  xtemp->Update(1.0, *x_, 0.0);
  epPermQTMatrix->Multiply(false, *xtemp, *x_);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Epetra_Map> LINALG::SOLVER::KrylovSolver::ExtractPermutationMap(
    const std::string solverSublist, const std::string mapName)
{
  // extract (user-given) additional information about linear system from
  // "Aztec/Belos Parameters" -> "Linear System properties"
  Teuchos::RCP<Epetra_Map> epSlaveDofMap = Teuchos::null;
  if (Params().isSublist(solverSublist) &&
      Params().sublist(solverSublist).isSublist("Linear System properties"))
  {
    Teuchos::ParameterList& solverParams = Params().sublist(solverSublist);
    Teuchos::ParameterList& linSystemProps = solverParams.sublist("Linear System properties");

    if (linSystemProps.isParameter(mapName))
      epSlaveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>(mapName);
  }
  return epSlaveDofMap;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
bool LINALG::SOLVER::KrylovSolver::DecideAboutPermutation(const Teuchos::RCP<Epetra_CrsMatrix>& A)
{
  bool bPermutationRecommended = false;

  // extract permuted matrix A (without scaling)
  Teuchos::RCP<Matrix> xPermA = data_->Get<Teuchos::RCP<Matrix>>("permA", PermFact_.get());

  // find problematic rows/columns in permuted and original matrix
  Teuchos::RCP<Map> PermutedNonDiagMap = FindNonDiagonalDominantRows(xPermA, diagDominanceRatio_);
  Teuchos::RCP<Map> NonPermutedNonDiagMap = FindNonDiagonalDominantRows(A, 1.0);

  double tolerance = 1e-5;
  Teuchos::RCP<Map> PermutedNearZeroMap = FindZeroDiagonalEntries(xPermA, tolerance);
  Teuchos::RCP<Map> NonPermutedNearZeroMap = FindZeroDiagonalEntries(A, tolerance);

  int PermutedNearZeros = PermutedNearZeroMap->getGlobalNumElements();
  int NonPermutedNearZeros = NonPermutedNearZeroMap->getGlobalNumElements();

  int NonPermutedZeros = CountZerosOnDiagonalEpetra(A);
  int PermutedZeros = CountZerosOnDiagonal(xPermA);

  Teuchos::RCP<Teuchos::FancyOStream> fos =
      Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(0);
  *fos << "---------------------------- MATRIX analysis -----------------------" << std::endl;
  *fos << "| permutation strategy: " << std::setw(25)
       << PermFact_->GetParameter("PermutationStrategy") << "                  |" << std::endl;
  *fos << "| non-permuted matrix A              | permuted matrix A           |" << std::endl;
  *fos << "| zeros on diagonal: " << std::setw(8) << std::right << NonPermutedZeros << "        ";
  *fos << "| zeros on diagonal: " << std::setw(8) << std::right << PermutedZeros << " |"
       << std::endl;
  *fos << "| #near zeros on diagonal: " << std::setw(8) << std::right
       << NonPermutedNearZeroMap->getGlobalNumElements() << "  ";
  *fos << "| #near zeros on diag: " << std::setw(6) << std::right
       << PermutedNearZeroMap->getGlobalNumElements() << " |" << std::endl;
  *fos << "| (tolerance: " << std::setw(8) << std::right << tolerance << ")              ";
  *fos << "| (tolerance: " << std::setw(8) << std::right << tolerance << ")       |" << std::endl;
  *fos << "| #nonDiagDomEntries  : " << std::setw(8) << std::right
       << NonPermutedNonDiagMap->getGlobalNumElements() << "     ";
  *fos << "| #nonDiagDomEntries  : " << std::setw(5) << std::right
       << PermutedNonDiagMap->getGlobalNumElements() << " |" << std::endl;
  if (NonPermutedZeros == 0 && PermutedZeros > 0)
  {
    *fos << "| permutation failed. use original matrix instead                  |" << std::endl;
    bPermutationRecommended = false;
  }
  else if (NonPermutedZeros > 0 && PermutedZeros > 0)
  {
    // TODO merge this if clause with the next if clause below...
    if (NonPermutedNonDiagMap->getGlobalNumElements() < PermutedNonDiagMap->getGlobalNumElements())
    {
      *fos << "| permutation failed.                                                |" << std::endl;
      *fos << "| original linear system seems to have less problematic rows.        |" << std::endl;
      *fos << "| -> use non-permuted original linear system                         |" << std::endl;
      bPermutationRecommended = false;
    }
    else
    {
      *fos << "| permutation failed.                                                |" << std::endl;
      *fos << "| permuted linear system seems to have less problematic rows.        |" << std::endl;
      *fos << "| -> use permuted original linear system                             |" << std::endl;
      bPermutationRecommended = true;
    }
  }
  else
  {
    if (NonPermutedNonDiagMap->getGlobalNumElements() < PermutedNonDiagMap->getGlobalNumElements())
    {
      *fos << "| original linear system seems to have less problematic rows.      |" << std::endl;
      *fos << "| -> use non-permuted original linear system                       |" << std::endl;
      bPermutationRecommended = false;
    }
    else
    {
      *fos << "| permuted linear system seems to have less problematic rows.      |" << std::endl;
      *fos << "| -> use permuted linear system                                    |" << std::endl;
      bPermutationRecommended = true;
    }
  }
  *fos << "---------------------------- MATRIX analysis -----------------------" << std::endl;

  // set data depending on decision whether the permuted or the original matrix
  // shall be used
  // used for output in permutedAztecSolver
  if (bPermutationRecommended)
  {
    data_->Set("nonDiagDomRows", Teuchos::as<int>(PermutedNonDiagMap->getGlobalNumElements()));
  }
  else
  {
    data_->Set("nonDiagDomRows", Teuchos::as<int>(NonPermutedNonDiagMap->getGlobalNumElements()));
  }
  data_->Set("NonPermutedZerosOnDiagonal", NonPermutedZeros);
  data_->Set("PermutedZerosOnDiagonal", PermutedZeros);
  data_->Set("PermutedNearZeros", PermutedNearZeros);
  data_->Set("NonPermutedNearZeros", NonPermutedNearZeros);

  // feed preconditioner with more information about linear system using
  // the "Linear System properties" sublist in the preconditioner's
  // paramter list
  if (Preconditioner() != NULL)
  {
    const std::string precondParamListName = Preconditioner()->getParameterListName();
    if (Params().isSublist(precondParamListName))
    {
      Teuchos::ParameterList& precondParams = Params().sublist(precondParamListName);
      Teuchos::ParameterList& linSystemProps = precondParams.sublist("Linear System properties");

      // set information which rows are not diagonal dominant in matrix
      // e.g. the MueLu_Contact2 preconditioner uses this information to mark
      // these problematic rows and use 1 point aggregates
      if (bPermutationRecommended)
      {
        linSystemProps.set<Teuchos::RCP<Map>>("non diagonal-dominant row map", PermutedNonDiagMap);
        linSystemProps.set<Teuchos::RCP<Map>>("near-zero diagonal row map", PermutedNearZeroMap);
      }
      else
      {
        linSystemProps.set<Teuchos::RCP<Map>>(
            "non diagonal-dominant row map", NonPermutedNonDiagMap);
        linSystemProps.set<Teuchos::RCP<Map>>("near-zero diagonal row map", NonPermutedNearZeroMap);
      }
    }
  }

  return bPermutationRecommended;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<const Epetra_CrsMatrix> LINALG::SOLVER::KrylovSolver::GetOperator(
    const std::string name, const Teuchos::RCP<FactoryBase>& fact)
{
  Teuchos::RCP<Matrix> xPermScalOp = data_->Get<Teuchos::RCP<Matrix>>(name, fact.get());
  Teuchos::RCP<CrsMatrixWrap> xPermScalCrsOp =
      Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xPermScalOp);
  Teuchos::RCP<CrsMatrix> xPermScalCrsMat = xPermScalCrsOp->getCrsMatrix();
  Teuchos::RCP<EpetraCrsMatrix> xEpPermScalCrsMat =
      Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(xPermScalCrsMat);
  return xEpPermScalCrsMat->getEpetra_CrsMatrix();
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Epetra_CrsMatrix> LINALG::SOLVER::KrylovSolver::GetOperatorNonConst(
    const std::string name, const Teuchos::RCP<FactoryBase>& fact)
{
  Teuchos::RCP<Matrix> xPermScalOp = data_->Get<Teuchos::RCP<Matrix>>(name, fact.get());
  Teuchos::RCP<CrsMatrixWrap> xPermScalCrsOp =
      Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xPermScalOp);
  Teuchos::RCP<CrsMatrix> xPermScalCrsMat = xPermScalCrsOp->getCrsMatrix();
  Teuchos::RCP<EpetraCrsMatrix> xEpPermScalCrsMat =
      Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(xPermScalCrsMat);
  return xEpPermScalCrsMat->getEpetra_CrsMatrixNonConst();
}
