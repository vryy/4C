/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of solver type base class

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_SOLVERTYPE_H_
#define SOLVER_SOLVERTYPE_H_

#include "lib_dserror.H"
#include "linalg_krylov_projector.H"

namespace LINALG::SOLVER
{
  /// linear solver type base class
  template <class MatrixType, class VectorType>
  class SolverType
  {
   public:
    //! Destructor
    virtual ~SolverType() = default;

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
        Teuchos::RCP<LINALG::KrylovProjector> projector) = 0;

    virtual int Solve() = 0;

    virtual int ApplyInverse(const VectorType& X, VectorType& Y) = 0;

    /// return number of iterations performed by solver
    virtual int getNumIters() const
    {
      dserror("Not implemented in base class!");
      return -1;
    };
  };
}  // namespace LINALG::SOLVER

#endif
