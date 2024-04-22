/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 0

*/
/*----------------------------------------------------------------------*/

#include <Amesos_Klu.h>
#include <Amesos_Lapack.h>
#ifdef HAVE_AMESOS_UMFPACK
#include <Amesos_Umfpack.h>
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include <Amesos_Superludist.h>
#endif

#include "baci_linalg_krylov_projector.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_linear_solver_method_direct.hpp"

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
CORE::LINEAR_SOLVER::DirectSolver<MatrixType, VectorType>::DirectSolver(std::string solvertype)
    : solvertype_(solvertype),
      factored_(false),
      x_(Teuchos::null),
      b_(Teuchos::null),
      A_(Teuchos::null),
      amesos_(Teuchos::null),
      reindexer_(Teuchos::null),
      projector_(Teuchos::null)
{
  lp_ = Teuchos::rcp(new Epetra_LinearProblem());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::DirectSolver<MatrixType, VectorType>::Setup(
    Teuchos::RCP<MatrixType> matrix, Teuchos::RCP<VectorType> x, Teuchos::RCP<VectorType> b,
    const bool refactor, const bool reset, Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
{
  // Assume the input matrix to be a single block matrix
  bool bIsCrsMatrix = true;

  // try to cast input matrix to a Epetra_CrsMatrix
  // if the cast fails, the input matrix is a blocked operator which cannot be handled by a direct
  // solver unless the matrix is merged. This is a very expensive operation
  Teuchos::RCP<Epetra_CrsMatrix> crsA = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(matrix);
  if (crsA == Teuchos::null) bIsCrsMatrix = false;

  // store projector in internal member variable
  // If no projector is used it is Teuchos::null
  projector_ = projector;

  // Set internal member variables (store system matrix, rhs vector and solution vector)
  if (projector_ != Teuchos::null)
  {
    // instead of
    //
    // A x = b
    //
    // solve
    //
    // P^T A P x_tilda = P^T b
    //

    // cast system matrix to CORE::LINALG::SparseMatrix
    // check whether cast was successfull
    if (crsA == Teuchos::null)
    {
      FOUR_C_THROW("Could not cast system matrix to Epetra_CrsMatrix.");
    }
    // get view on systemmatrix as CORE::LINALG::SparseMatrix - this is no copy!
    CORE::LINALG::SparseMatrix A_view(crsA, CORE::LINALG::View);

    // apply projection to A without computing projection matrix thus avoiding
    // matrix-matrix multiplication
    Teuchos::RCP<CORE::LINALG::SparseMatrix> A2 = projector_->Project(A_view);

    // hand matrix over to A_
    A_ = A2->EpetraMatrix();
    // hand over to b_ and project to (P^T b)
    b_ = b;
    projector_->ApplyPT(*b_);
    // hand over x_ as x_tilda (only zeros yet, )
    x_ = x;
  }
  else
  {
    if (bIsCrsMatrix == true)
    {
      A_ = matrix;
    }
    else
    {
      Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> Ablock =
          Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(matrix);

      int matrixDim = Ablock->FullRangeMap().NumGlobalElements();
      if (matrixDim > 50000)
      {
        Teuchos::RCP<Teuchos::FancyOStream> fos =
            Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
        fos->setOutputToRootOnly(0);
        *fos << "---------------------------- ATTENTION -----------------------" << std::endl;
        *fos << "  Merging a " << Ablock->Rows() << " x " << Ablock->Cols() << " block matrix "
             << std::endl;
        *fos << "  of size " << matrixDim << " is a very expensive operation. " << std::endl;
        *fos << "  For performance reasons, try iterative solvers instead!" << std::endl;
        *fos << "---------------------------- ATTENTION -----------------------" << std::endl;
      }

      Teuchos::RCP<CORE::LINALG::SparseMatrix> Ablock_merged = Ablock->Merge();
      A_ = Ablock_merged->EpetraMatrix();
    }
    x_ = x;
    b_ = b;
  }

  // fill the linear problem
  lp_->SetRHS(b_.get());
  lp_->SetLHS(x_.get());
  lp_->SetOperator(A_.get());

  /* update reindexing of vectors: Doing so, we don't need to reset the solver
   * to enforce reindexing of the entire Epetra_LinearProblem. This allows for
   * reuse of the factorization.
   */
  if (not reindexer_.is_null() and not(reset or refactor)) reindexer_->fwd();

  if (reset or refactor or not IsFactored())
  {
    amesos_ = Teuchos::null;

    reindexer_ = Teuchos::rcp(new EpetraExt::LinearProblem_Reindex2(nullptr));

    if (solvertype_ == "umfpack")
    {
#ifdef HAVE_AMESOS_UMFPACK
      amesos_ = Teuchos::rcp(new Amesos_Umfpack((*reindexer_)(*lp_)));
#else
      FOUR_C_THROW(
          "UMFPACK was chosen as linear solver, but is not available in the configuration!");
#endif
    }
    else if (solvertype_ == "superlu")
    {
#ifdef HAVE_AMESOS_SUPERLUDIST
      amesos_ = Teuchos::rcp(new Amesos_Superludist((*reindexer_)(*lp_)));
#else
      FOUR_C_THROW(
          "Superlu was chosen as linear solver, but is not available in the configuration!");
#endif
    }
    else
    {
      FOUR_C_THROW("UMFPACK or Superlu have to be available to use a direct linear solver.");
    }

    factored_ = false;
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
int CORE::LINEAR_SOLVER::DirectSolver<MatrixType, VectorType>::Solve()
{
  if (amesos_ == Teuchos::null) FOUR_C_THROW("No solver allocated");

  // Problem has not been factorized before
  if (not IsFactored())
  {
    int err = amesos_->SymbolicFactorization();
    if (err) FOUR_C_THROW("Amesos::SymbolicFactorization returned an err");
    err = amesos_->NumericFactorization();
    if (err) FOUR_C_THROW("Amesos::NumericFactorization returned an err");

    factored_ = true;
  }

  int err = amesos_->Solve();

  if (err) FOUR_C_THROW("Amesos::Solve returned an err");

  if (projector_ != Teuchos::null)
  {
    // get x from x = P x_tilda
    projector_->ApplyP(*x_);
  }
  // direct solver does not support errorcodes
  return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// explicit initialization
template class CORE::LINEAR_SOLVER::DirectSolver<Epetra_Operator, Epetra_MultiVector>;

FOUR_C_NAMESPACE_CLOSE
