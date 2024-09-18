/*----------------------------------------------------------------------*/
/*! \file

\brief Description

\level 0

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_METHOD_DIRECT_HPP
#define FOUR_C_LINEAR_SOLVER_METHOD_DIRECT_HPP

#include "4C_config.hpp"

#include "4C_linear_solver_method.hpp"

#include <Amesos_BaseSolver.h>
#include <EpetraExt_Reindex_LinearProblem2.h>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinearSolver
{
  /// direct linear solver (using amesos)
  template <class MatrixType, class VectorType>
  class DirectSolver : public SolverTypeBase<MatrixType, VectorType>
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
    void setup(Teuchos::RCP<MatrixType> matrix, Teuchos::RCP<VectorType> x,
        Teuchos::RCP<VectorType> b, const bool refactor, const bool reset,
        Teuchos::RCP<Core::LinAlg::KrylovProjector> projector = Teuchos::null) override;

    //! Actual call to the underlying amesos solver
    int solve() override;

    bool is_factored() { return factored_; }

   private:
    //! type/implementation of Amesos solver to be used
    const std::string solvertype_;

    //! flag indicating whether a valid factorization is stored
    bool factored_;

    //! a linear problem wrapper class used by Trilinos and for scaling of the system
    Teuchos::RCP<Epetra_LinearProblem> linear_problem_;

    //! initial guess and solution
    Teuchos::RCP<VectorType> x_;

    //! right hand side vector
    Teuchos::RCP<VectorType> b_;

    //! system of equations
    Teuchos::RCP<MatrixType> a_;

    //! an abstract amesos solver that can be any of the amesos concrete implementations
    Teuchos::RCP<Amesos_BaseSolver> solver_;

    //! reindex linear problem for amesos
    Teuchos::RCP<EpetraExt::LinearProblem_Reindex2> reindexer_;

    /*! \brief Krylov projector for solving near singular linear systems
     *
     * Instead of solving Ax=b a projected system of the form P'APu=P'b is solved.
     * With P being the krylov projector.
     *
     * P. Bochev and R. B. Lehoucq: On the Finite Element Solution of the Pure Neumann Problem,
     * SIAM Review, 47(1):50-66, 2005, http://dx.doi.org/10.1137/S0036144503426074
     *
     */
    Teuchos::RCP<Core::LinAlg::KrylovProjector> projector_;
  };
}  // namespace Core::LinearSolver

FOUR_C_NAMESPACE_CLOSE

#endif
