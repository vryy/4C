/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of interface to Belos solver package

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_BELOSSOLVER_H_
#define SOLVER_BELOSSOLVER_H_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "solver_krylovsolver.H"
#include "solver_preconditionertype.H"

namespace LINALG::SOLVER
{
  /// krylov subspace linear solvers (belos) with right-side preconditioning
  template <class MatrixType, class VectorType>
  class BelosSolver : public KrylovSolver<MatrixType, VectorType>
  {
   public:
    //! Constructor
    BelosSolver(const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile);

    /*! \brief Setup the solver object
     *
     * @param A Matrix of the linear system
     * @param x Solution vector of the linear system
     * @param b Right-hand side vector of the linear system
     * @param refactor Boolean flag to enforce a refactorization of the matrix
     * @param reset Boolean flag to enforce a full reset of the solver object
     * @param projector Krylov projector
     */
    virtual void Setup(Teuchos::RCP<MatrixType> A, Teuchos::RCP<VectorType> x,
        Teuchos::RCP<VectorType> b, const bool refactor, const bool reset,
        Teuchos::RCP<LINALG::KrylovProjector> projector);

    //! Actual call to the underlying Belos solver
    virtual int Solve();
  };
}  // namespace LINALG::SOLVER

#endif /* SOLVER_BELOSSOLVER_H_ */
