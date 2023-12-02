/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of interface to Belos solver package

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_METHOD_BELOS_H
#define BACI_LINEAR_SOLVER_METHOD_BELOS_H

#include "baci_linear_solver_method_krylov.H"
#include "baci_linear_solver_preconditioner_type.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace CORE::LINEAR_SOLVER
{
  /// krylov subspace linear solvers (belos) with right-side preconditioning
  template <class MatrixType, class VectorType>
  class BelosSolver : public KrylovSolver<MatrixType, VectorType>
  {
   public:
    //! Constructor
    BelosSolver(const Epetra_Comm& comm, Teuchos::ParameterList& params);

    /*! \brief Setup the solver object
     *
     * @param A Matrix of the linear system
     * @param x Solution vector of the linear system
     * @param b Right-hand side vector of the linear system
     * @param refactor Boolean flag to enforce a refactorization of the matrix
     * @param reset Boolean flag to enforce a full reset of the solver object
     * @param projector Krylov projector
     */
    void Setup(Teuchos::RCP<MatrixType> matrix, Teuchos::RCP<VectorType> x,
        Teuchos::RCP<VectorType> b, const bool refactor, const bool reset,
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector) override;

    //! Actual call to the underlying Belos solver
    int Solve() override;

    //! return number of iterations
    int getNumIters() const override { return numiters_; };

   private:
    //! number of iterations
    int numiters_;
  };
}  // namespace CORE::LINEAR_SOLVER

#endif
