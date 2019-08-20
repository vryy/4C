/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of interface to Belos solver package

\level 1

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_BELOSSOLVER_H_
#define SOLVER_BELOSSOLVER_H_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "solver_krylovsolver.H"
#include "solver_preconditionertype.H"

namespace LINALG
{
  class KrylovProjector;

  namespace SOLVER
  {
    /// krylov subspace linear solvers (belos) with right-side preconditioning
    class BelosSolver : public KrylovSolver
    {
     public:
      //! Constructor
      BelosSolver(const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile);

      //! Destructor
      virtual ~BelosSolver();

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
          Teuchos::RCP<LINALG::KrylovProjector> projector);

      //! Actual call to the underlying Belos solver
      virtual int Solve();
    };
  }  // namespace SOLVER
}  // namespace LINALG


#endif /* SOLVER_BELOSSOLVER_H_ */
