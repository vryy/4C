/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration

\brief Declaration
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/

// Amesos headers
#include <Amesos_Klu.h>
#include <Amesos_Lapack.h>
#ifdef HAVE_AMESOS_UMFPACK
#include <Amesos_Umfpack.h>
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include <Amesos_Superludist.h>
#endif

// EpetraExt headers
#include <EpetraExt_Reindex_LinearProblem2.h>

// BACI headers
#include "solver_directsolver.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_krylov_projector.H"
#include <Epetra_CrsMatrix.h>
#include <Teuchos_Time.hpp>

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::DirectSolver::DirectSolver(std::string solvertype)
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
LINALG::SOLVER::DirectSolver::~DirectSolver()
{
  amesos_ = Teuchos::null;
  reindexer_ = Teuchos::null;
  lp_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::DirectSolver::Setup(Teuchos::RCP<Epetra_Operator> matrix,
    Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> b, const bool refactor,
    const bool reset, Teuchos::RCP<LINALG::KrylovProjector> projector)
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

    // cast system matrix to LINALG::SparseMatrix
    // check whether cast was successfull
    if (crsA == Teuchos::null)
    {
      dserror("Could not cast system matrix to Epetra_CrsMatrix.");
    }
    // get view on systemmatrix as LINALG::SparseMatrix - this is no copy!
    LINALG::SparseMatrix A_view(crsA, View);

    // apply projection to A without computing projection matrix thus avoiding
    // matrix-matrix multiplication
    Teuchos::RCP<LINALG::SparseMatrix> A2 = projector_->Project(A_view);

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
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> Ablock =
          Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(matrix);

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

      Teuchos::RCP<LINALG::SparseMatrix> Ablock_merged = Ablock->Merge();
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

    reindexer_ = Teuchos::rcp(new EpetraExt::LinearProblem_Reindex2(NULL));

    if (solvertype_ == "klu")
    {
      amesos_ = Teuchos::rcp(new Amesos_Klu((*reindexer_)(*lp_)));
    }
    else if (solvertype_ == "umfpack")
    {
#ifdef HAVE_AMESOS_UMFPACK
      amesos_ = Teuchos::rcp(new Amesos_Umfpack((*reindexer_)(*lp_)));
#else
      dserror("no umfpack here");
#endif
    }
    else if (solvertype_ == "superlu")
    {
#ifdef HAVE_AMESOS_SUPERLUDIST
      amesos_ = Teuchos::rcp(new Amesos_Superludist((*reindexer_)(*lp_)));
#else
      Teuchos::RCP<Teuchos::FancyOStream> fos =
          Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
      fos->setOutputToRootOnly(0);
      *fos << "Warning: No SuperLU_dist available. Linear system will be solved with KLU instead..."
           << std::endl;
      amesos_ = Teuchos::rcp(new Amesos_Klu((*reindexer_)(*lp_)));
#endif
    }
    else if (solvertype_ == "lapack")
    {
      amesos_ = Teuchos::rcp(new Amesos_Lapack((*reindexer_)(*lp_)));
    }

    factored_ = false;
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::DirectSolver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)
{
  x_->Update(1., X, 0.);
  Solve();
  Y.Update(1., *b_, 0.);

  return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::DirectSolver::Solve()
{
  if (amesos_ == Teuchos::null) dserror("No solver allocated");

  // Problem has not been factorized before
  if (not IsFactored())
  {
    int err = amesos_->SymbolicFactorization();
    if (err) dserror("Amesos::SymbolicFactorization returned an err");
    err = amesos_->NumericFactorization();
    if (err) dserror("Amesos::NumericFactorization returned an err");

    factored_ = true;
  }

  int err = amesos_->Solve();

  if (err) dserror("Amesos::Solve returned an err");

  if (projector_ != Teuchos::null)
  {
    // get x from x = P x_tilda
    projector_->ApplyP(*x_);
  }
  // direct solver does not support errorcodes
  return 0;
}
