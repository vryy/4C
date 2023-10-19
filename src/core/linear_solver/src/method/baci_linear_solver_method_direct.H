/*----------------------------------------------------------------------*/
/*! \file

\brief Description

\level 0

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_METHOD_DIRECT_H
#define BACI_LINEAR_SOLVER_METHOD_DIRECT_H

#include "baci_linear_solver_method.H"

#include <Amesos_BaseSolver.h>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_Reindex_LinearProblem2.h>

namespace CORE::LINEAR_SOLVER
{
  /// direct linear solver (using amesos)
  template <class MatrixType, class VectorType>
  class DirectSolver : public SolverType<MatrixType, VectorType>
  {
   public:
    //! Constructor
    explicit DirectSolver(std::string solvertype);

    /*! \brief Setup the solver object
     *
     * @param A Matrix of the linear system
     * @param x Solution vector of the linear system
     * @param b Right-hand side vector of the linear system
     * @param refactor Boolean flag to enforce a refactorization of the matrix
     * @param reset Boolean flag to enforce a full reset of the solver object
     * @param projector Krylov projector
     */
    void Setup(Teuchos::RCP<MatrixType> A, Teuchos::RCP<VectorType> x, Teuchos::RCP<VectorType> b,
        const bool refactor, const bool reset,
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector = Teuchos::null) override;

    //! Actual call to the underlying amesos solver
    int Solve() override;

    int ApplyInverse(const VectorType& X, VectorType& Y) override;

    bool IsFactored() { return factored_; }

   private:
    //! type/implementation of Amesos solver to be used
    const std::string solvertype_;

    //! flag indicating whether a valid factorization is stored
    bool factored_;

    //! a linear problem wrapper class used by Trilinos and for scaling of the system
    Teuchos::RCP<Epetra_LinearProblem> lp_;

    //! initial guess and solution
    Teuchos::RCP<VectorType> x_;

    //! right hand side vector
    Teuchos::RCP<VectorType> b_;

    //! system of equations
    Teuchos::RCP<MatrixType> A_;

    //! an abstract amesos solver that can be any of the amesos concrete implementations
    Teuchos::RCP<Amesos_BaseSolver> amesos_;

    //! reindex linear problem for amesos
    Teuchos::RCP<EpetraExt::LinearProblem_Reindex2> reindexer_;

    //! Krylov projector if necessary
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector_;
  };
}  // namespace CORE::LINEAR_SOLVER

#endif
