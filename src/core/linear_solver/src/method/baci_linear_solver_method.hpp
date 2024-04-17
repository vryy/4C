/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of solver type base class

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_METHOD_HPP
#define FOUR_C_LINEAR_SOLVER_METHOD_HPP

#include "baci_config.hpp"

#include "baci_linalg_krylov_projector.hpp"
#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
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
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector) = 0;

    virtual int Solve() = 0;

    /// return number of iterations performed by solver
    virtual int getNumIters() const
    {
      dserror("Not implemented in base class!");
      return -1;
    };
  };
}  // namespace CORE::LINEAR_SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
