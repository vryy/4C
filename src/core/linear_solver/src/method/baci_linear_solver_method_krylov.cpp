/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of Baci's interface to Krylov solvers

\level 0

*/
/*---------------------------------------------------------------------*/

#include "baci_linear_solver_method_krylov.H"

#include "baci_linear_solver_amgnxn_preconditioner.H"
#include "baci_linear_solver_preconditioner_block.H"
#include "baci_linear_solver_preconditioner_ifpack.H"
#include "baci_linear_solver_preconditioner_krylovprojection.H"
#include "baci_linear_solver_preconditioner_ml.H"
#include "baci_linear_solver_preconditioner_muelu.H"
#include "baci_linear_solver_preconditioner_point.H"
#include "baci_utils_exceptions.H"

#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <MueLu.hpp>
#include <MueLu_FactoryBase.hpp>
#include <MueLu_HierarchyUtils.hpp>
#include <MueLu_PermutationFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_SmootherPrototype.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_VerboseObject.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

using Map = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
using CrsMatrix = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
using CrsMatrixWrap = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
using Vector = Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
using MultiVector = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
using VectorFactory = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
using MapFactory = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>;

using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

using FactoryManager = MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
using FactoryBase = MueLu::FactoryBase;

using SC = Scalar;
using LO = LocalOrdinal;
using GO = GlobalOrdinal;
using NO = Node;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::KrylovSolver(
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
template <class MatrixType, class VectorType>
bool CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::AllowReusePreconditioner(
    const int reuse, const bool reset)
{
  bool bAllowReuse = true;  // default: allow reuse of preconditioner

  // first, check some parameters with information that has to be updated
  const Teuchos::ParameterList& linSysParams = Params().sublist("Belos Parameters");

  CheckReuseStatusOfActiveSet(bAllowReuse, &linSysParams);

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
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::CheckReuseStatusOfActiveSet(
    bool& bAllowReuse, const Teuchos::ParameterList* linSysParams)
{
  if (linSysParams != nullptr)
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
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::CreatePreconditioner(
    Teuchos::ParameterList& solverlist, const bool isCrsMatrix,
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::Solver:  1.1)   CreatePreconditioner");

  preconditioner_ = Teuchos::null;

  if (isCrsMatrix)
  {
    // get type of preconditioner and build either Ifpack or ML
    // if we have an ifpack parameter list, we do ifpack
    // if we have an ml parameter list we do ml
    if (Params().isSublist("IFPACK Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::IFPACKPreconditioner(
          outfile_, Params().sublist("IFPACK Parameters"), solverlist));
    }
    else if (Params().isSublist("ML Parameters"))
    {
      preconditioner_ = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::MLPreconditioner(outfile_, Params().sublist("ML Parameters")));
    }
    else if (Params().isSublist("MueLu Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuPreconditioner(
          outfile_, Params().sublist("MueLu Parameters")));
    }
    else if (Params().isSublist("MueLu (BeamSolid) Parameters"))
    {
#ifdef BACI_WITH_TRILINOS_DEVELOP
      preconditioner_ = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::MueLuBeamSolidBlockPreconditioner(outfile_, Params()));
#else
      dserror("MueLu (BeamSolid) preconditioner only available in Trilinos_Develop.");
#endif
    }
    else
    {
      dserror("unknown preconditioner");
    }

    // decide whether we do what kind of scaling
    std::string scaling = solverlist.get("scaling", "none");
    if (scaling == "none")
    {
    }
    else if (scaling == "infnorm")
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::InfNormPreconditioner(preconditioner_));
    }
    else if (scaling == "symmetric")
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::SymDiagPreconditioner(preconditioner_));
    }
    else
      dserror("Unknown type of scaling found in parameter list");

    if (projector != Teuchos::null)
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner(
          outfile_, preconditioner_, projector));
    }
  }
  else
  {
    // assume block matrix

    if (Params().isSublist(
            "SIMPLER"))  // old BACI::(Cheap)SIMPLER preconditioner TODO: remove/replace me
    {
      dserror("SIMPLER sublist not supported any more.");
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::SimplePreconditioner(
          outfile_, Params()));  // Michael's SIMPLE for Fluid
    }
    else if (Params().isSublist("CheapSIMPLE Parameters"))
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::SimplePreconditioner(outfile_, Params()));
    }
    else if (Params().isSublist("BGS Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::BGSPreconditioner(
          outfile_, Params(), Params().sublist("BGS Parameters")));
    }
    else if (Params().isSublist("MueLu (Fluid) Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuFluidBlockPreconditioner(
          outfile_, Params().sublist("MueLu (Fluid) Parameters")));
    }
    else if (Params().isSublist("MueLu (TSI) Parameters"))
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuTsiBlockPreconditioner(outfile_, Params()));
    }
    else if (Params().isSublist("MueLu (Contact) Parameters"))
    {
      preconditioner_ = Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuContactSpPreconditioner(
          outfile_, Params().sublist("MueLu (Contact) Parameters")));
    }
    else if (Params().isSublist("MueLu (FSI) Parameters"))
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuFsiBlockPreconditioner(outfile_, Params()));
    }
    else if (Params().isSublist("AMGnxn Parameters"))
    {
      preconditioner_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::AMGnxn_Preconditioner(outfile_, Params()));
    }
    else
    {
#ifdef DEBUG
      Params().print(std::cout);
#endif

      dserror("unknown preconditioner for block matrix solver");
    }
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::BuildPermutationOperator(
    const Teuchos::RCP<Epetra_CrsMatrix>& A, const Teuchos::RCP<Epetra_Map>& epSlaveDofMap)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new EpetraCrsMatrix(A));
  Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>> xCrsOp =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>(xCrsA));
  Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>> xOp =
      Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LO, GO, Node>>(xCrsOp);
  xOp->SetFixedBlockSize(this->Params()
                             .sublist("NodalBlockInformation")
                             .template get<int>("number of momentum dofs"));  // set nBlockSize

  data_->setDefaultVerbLevel(Teuchos::VERB_NONE);
  data_->setlib(Xpetra::UseEpetra);
  data_->Set("A", xOp);


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // check, if "SlaveDofMap" information is available in parameter lists
  if (epSlaveDofMap != Teuchos::null)
  {
    Teuchos::RCP<Xpetra::EpetraMapT<GO, NO>> xSlaveDofMap =
        Teuchos::rcp(new Xpetra::EpetraMapT<GO, NO>(epSlaveDofMap));

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
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::PermuteLinearSystem(
    const Teuchos::RCP<Epetra_CrsMatrix>& A, const Teuchos::RCP<Epetra_MultiVector>& b)
{
  if (!data_->IsAvailable("A", PermFact_.get()) || !data_->IsAvailable("permP", PermFact_.get()) ||
      !data_->IsAvailable("permScaling", PermFact_.get()))
    dserror("PermutedIterativeSolver: call BuildPermutationOperator before PermuteLinearSystem");

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
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::PermuteNullSpace(
    const Teuchos::RCP<Epetra_CrsMatrix>& A)
{
  // note: we usually do not permute the null space!!!

  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xOp = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xOp->SetFixedBlockSize(this->Params()
                             .sublist("NodalBlockInformation")
                             .template get<int>("number of momentum dofs"));  // set nBlockSize

  // detect MueLu/ML Paramter list
  std::string MultiGridParameterListName = "";
  if (Params().isSublist("ML Parameters"))
    MultiGridParameterListName = "ML Parameters";
  else if (Params().isSublist("MueLu Parameters"))
    MultiGridParameterListName = "MueLu Parameters";
  else if (Params().isSublist("MueLu (Contact) Parameters"))
    MultiGridParameterListName = "MueLu (Contact) Parameters";

  // retransform nullspace vectors
  if (MultiGridParameterListName == "") return;  // no nullspace to permute

  Teuchos::RCP<Matrix> xPermQtMatrix = data_->Get<Teuchos::RCP<Matrix>>("permQT", PermFact_.get());
  int numdf =
      this->Params().sublist(MultiGridParameterListName).template get<int>("PDE equations", -1);
  int dimns = this->Params()
                  .sublist(MultiGridParameterListName)
                  .template get<int>("null space: dimension", -1);
  if (dimns == -1 || numdf == -1)
    dserror(
        "PermutedIterativeSolver: Error in MueLu/ML parameters: PDE equations or null space "
        "dimension "
        "wrong.");
  Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> rowMap = xOp->getRowMap();

  Teuchos::RCP<MultiVector> nspVector =
      Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
          rowMap, dimns, true);
  Teuchos::RCP<std::vector<double>> nsdata =
      this->Params()
          .sublist(MultiGridParameterListName)
          .template get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);

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
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
Teuchos::RCP<Map>
CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::FindNonDiagonalDominantRows(
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
  const int numLocalElements = xA->getRowMap()->getLocalNumElements();

  for (int row = 0; row < numLocalElements; row++)
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
      NonDiagonalDominantGIDs.data(), NonDiagonalDominantGIDs.size());

  Teuchos::RCP<Map> NonDiagonalDominantGIDsMap = MapFactory::Build(xA->getRowMap()->lib(),
      Teuchos::OrdinalTraits<int>::invalid(), NonDiagonalDominantGIDs_view, 0, comm);

  return NonDiagonalDominantGIDsMap;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
Teuchos::RCP<Map>
CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::FindNonDiagonalDominantRows(
    const Teuchos::RCP<Epetra_CrsMatrix>& A, double diagDominanceRatio)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xA = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xA->SetFixedBlockSize(this->Params()
                            .sublist("NodalBlockInformation")
                            .template get<int>("number of momentum dofs"));  // set nBlockSize

  return FindNonDiagonalDominantRows(xA, diagDominanceRatio);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
Teuchos::RCP<Map>
CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::FindZeroDiagonalEntries(
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

  const int numLocalElements = xA->getRowMap()->getLocalNumElements();

  for (int i = 0; i < numLocalElements; ++i)
  {
    if (std::abs(diagAVecData[i]) < tolerance)
    {  // pick out all rows with very small diagonal entries
      lNumZeros++;
      zeroGids.push_back(diagAVec->getMap()->getGlobalElement(i));
    }
  }

  const Teuchos::ArrayView<const LocalOrdinal> zeroGids_view(zeroGids.data(), zeroGids.size());

  Teuchos::RCP<Map> zeroDiagonalMap = MapFactory::Build(
      xA->getRowMap()->lib(), Teuchos::OrdinalTraits<int>::invalid(), zeroGids_view, 0, comm);

  return zeroDiagonalMap;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
Teuchos::RCP<Map>
CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::FindZeroDiagonalEntries(
    const Teuchos::RCP<Epetra_CrsMatrix>& A, double tolerance)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xA = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xA->SetFixedBlockSize(this->Params()
                            .sublist("NodalBlockInformation")
                            .template get<int>("number of momentum dofs"));  // set nBlockSize

  return FindZeroDiagonalEntries(xA, tolerance);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
int CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::CountZerosOnDiagonalEpetra(
    const Teuchos::RCP<Epetra_CrsMatrix>& A)
{
  // wrap Epetra_CrsMatrix -> Xpetra::Matrix
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> xCrsA =
      Teuchos::rcp(new EpetraCrsMatrix(A));
  Teuchos::RCP<CrsMatrixWrap> xCrsOp = Teuchos::rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<Matrix> xOp = Teuchos::rcp_dynamic_cast<Matrix>(xCrsOp);
  xOp->SetFixedBlockSize(this->Params()
                             .sublist("NodalBlockInformation")
                             .template get<int>("number of momentum dofs"));  // set nBlockSize

  return CountZerosOnDiagonal(xOp);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
int CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::CountZerosOnDiagonal(
    const Teuchos::RCP<const Matrix>& xOp)
{
  Teuchos::RCP<Vector> diagAVec = VectorFactory::Build(xOp->getRowMap(), true);
  xOp->getLocalDiagCopy(*diagAVec);
  Teuchos::ArrayRCP<const Scalar> diagAVecData = diagAVec->getData(0);
  LocalOrdinal lNumZeros = 0;
  GlobalOrdinal gNumZeros = 0;

  const int numLocalElements = diagAVec->getMap()->getLocalNumElements();

  for (int i = 0; i < numLocalElements; ++i)
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
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::ReTransformSolution()
{
  Teuchos::RCP<const Epetra_CrsMatrix> epPermQTMatrix = GetOperator("permQT", PermFact_);
  Teuchos::RCP<Epetra_MultiVector> xtemp = Teuchos::rcp(new Epetra_MultiVector(*x_));
  xtemp->Update(1.0, *x_, 0.0);
  epPermQTMatrix->Multiply(false, *xtemp, *x_);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
Teuchos::RCP<Epetra_Map>
CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::ExtractPermutationMap(
    const std::string solverSublist, const std::string mapName)
{
  // extract (user-given) additional information about linear system from
  // "Belos Parameters" -> "Linear System properties"
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
template <class MatrixType, class VectorType>
bool CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::DecideAboutPermutation(
    const Teuchos::RCP<Epetra_CrsMatrix>& A)
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
  // used for output in permutedIterativeSolver
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
  {
    const std::string precondParamListName = Preconditioner().getParameterListName();
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
template <class MatrixType, class VectorType>
Teuchos::RCP<const Epetra_CrsMatrix>
CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::GetOperator(
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
template <class MatrixType, class VectorType>
Teuchos::RCP<Epetra_CrsMatrix>
CORE::LINEAR_SOLVER::KrylovSolver<MatrixType, VectorType>::GetOperatorNonConst(
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

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// explicit initialization
template class CORE::LINEAR_SOLVER::KrylovSolver<Epetra_Operator, Epetra_MultiVector>;
