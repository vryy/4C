/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of interface to AztecOO solver package

\level 0

*/
/*---------------------------------------------------------------------*/
#ifndef SOLVER_AZTECSOLVER_H_
#define SOLVER_AZTECSOLVER_H_

#include <MueLu_PermutationFactory_fwd.hpp>
#include <MueLu_FactoryBase_fwd.hpp>
#include <MueLu_Level_fwd.hpp>
#include <MueLu_UseDefaultTypes.hpp>

#include <AztecOO_StatusTest.h>
#include <AztecOO_StatusTestCombo.h>
#include <AztecOO_StatusTestMaxIters.h>
#include <AztecOO_StatusTestResNorm.h>

#include "solver_aztecsolver_projectedresidual.H"
#include "solver_krylovsolver.H"

namespace LINALG::SOLVER
{
  /// krylov subspace linear solvers (aztec) with right-side preconditioning
  template <class MatrixType, class VectorType>
  class AztecSolver : public KrylovSolver<MatrixType, VectorType>
  {
   public:
    //! Constructor
    AztecSolver(const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile);

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
        Teuchos::RCP<LINALG::KrylovProjector> projector) override;

    //! Actual call to the underlying Aztec solver
    int Solve() override;

    //! return number of iterations
    int getNumIters() const override { return numiters_; };

   private:
    //! the Krylov projector
    Teuchos::RCP<LINALG::KrylovProjector> projector_;
    //! an aztec maxiter status test
    Teuchos::RCP<AztecOO_StatusTestMaxIters> aztest_maxiter_;
    //! an aztec projected residual 2-norm convergence test
    Teuchos::RCP<AztecOO_StatusTestProjResNorm> aztest_norm2_;
    //! a test combination (for norm2 and maxiter)
    Teuchos::RCP<AztecOO_StatusTestCombo> aztest_combo_;
    //! number of iterations
    int numiters_;
  };
}  // namespace LINALG::SOLVER

#endif
