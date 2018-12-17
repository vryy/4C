/*!----------------------------------------------------------------------
\file solver_directsolver.H

<pre>
\brief Description
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_DIRECTSOLVER_H_
#define SOLVER_DIRECTSOLVER_H_

#include "solver_solvertype.H"

// forward declarations
class Amesos_BaseSolver;
namespace EpetraExt
{
  class LinearProblem_Reindex2;
}

namespace LINALG
{
  class KrylovProjector;

  namespace SOLVER
  {
    /// direct linear solver (using amesos)
    class DirectSolver : public SolverType
    {
     public:
      //! Constructor
      explicit DirectSolver(std::string solvertype);

      //! Desctructor
      virtual ~DirectSolver();

      //! Setup
      virtual void Setup(Teuchos::RCP<Epetra_Operator> matrix,  ///< matrix
          Teuchos::RCP<Epetra_MultiVector> x,                   ///< solution vector
          Teuchos::RCP<Epetra_MultiVector> b,                   ///< right hand side
          const bool refactor,                                  ///< Enforce refactorization?
          const bool reset,                                     ///< Enforce re-setup of solver?
          Teuchos::RCP<LINALG::KrylovProjector> projector = Teuchos::null  ///< Krylov projector
      );

      //! Actual call to the underlying amesos solver
      virtual int Solve();

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y);

      virtual PreconditionerType* Preconditioner() { return NULL; }

      bool IsFactored() { return factored_; }

     private:
      //! type/implementation of Amesos solver to be used
      const std::string solvertype_;

      //! flag indicating whether a valid factorization is stored
      bool factored_;

      //! a linear problem wrapper class used by Trilinos and for scaling of the system
      Teuchos::RCP<Epetra_LinearProblem> lp_;

      //! initial guess and solution
      Teuchos::RCP<Epetra_MultiVector> x_;

      //! right hand side vector
      Teuchos::RCP<Epetra_MultiVector> b_;

      //! system of equations
      Teuchos::RCP<Epetra_Operator> A_;

      //! an abstract amesos solver that can be any of the amesos concrete implementations
      Teuchos::RCP<Amesos_BaseSolver> amesos_;

      //! reindex linear problem for amesos
      Teuchos::RCP<EpetraExt::LinearProblem_Reindex2> reindexer_;

      //! Krylov projector if necessary
      Teuchos::RCP<LINALG::KrylovProjector> projector_;
    };
  }  // namespace SOLVER
}  // namespace LINALG


#endif /* SOLVER_DIRECTSOLVER_H_ */
