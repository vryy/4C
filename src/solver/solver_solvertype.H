/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of solver type base class

\level 1

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_SOLVERTYPE_H_
#define SOLVER_SOLVERTYPE_H_

#include "../drt_lib/drt_dserror.H"

namespace LINALG
{
  class KrylovProjector;
  namespace SOLVER
  {
    class PreconditionerType;

    /// linear solver type base class
    class SolverType
    {
     public:
      //! Destructor
      virtual ~SolverType() {}

      /*! \brief Setup the solver object
       *
       * @param matrix Matrix of the linear system
       * @param x Solution vector of the linear system
       * @param b Right-hand side vector of the linear system
       * @param refactor Boolean flag to enforce a refactorization of the matrix
       * @param reset Boolean flag to enforce a full reset of the solver object
       * @param projector Krylov projector
       */
      virtual void Setup(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
          Teuchos::RCP<Epetra_MultiVector> b, const bool refactor, const bool reset,
          Teuchos::RCP<LINALG::KrylovProjector> projector) = 0;

      virtual int Solve() = 0;

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) = 0;

      /// preconditioner used within this solver
      virtual PreconditionerType* Preconditioner() = 0;

      /// return number of iterations performed by solver
      virtual int getNumIters() const
      {
        dserror("Not implemented in base class!");
        return -1;
      };
    };
  }  // namespace SOLVER
}  // namespace LINALG


#endif /* SOLVER_SOLVERTYPE_H_ */
