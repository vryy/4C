/*!----------------------------------------------------------------------
\file solver_directsolver.H

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_DIRECTSOLVER_H_
#define SOLVER_DIRECTSOLVER_H_

#include "solver_solvertype.H"

// forward declarations
class Amesos_BaseSolver;
namespace EpetraExt {
  class LinearProblem_Reindex;
}

namespace LINALG
{
  namespace SOLVER
  {
    /// direct linear solver (using amesos)
    class DirectSolver : public SolverType
    {
    public:

      explicit DirectSolver( std::string solvertype );

      virtual ~DirectSolver();

      virtual void Setup( Teuchos::RCP<Epetra_Operator>          matrix        ,
    		  	  	  	  Teuchos::RCP<Epetra_MultiVector>       x             ,
    		  	  	  	  Teuchos::RCP<Epetra_MultiVector>       b             ,
                          bool                                   refactor      ,
                          bool                                   reset         ,
                          Teuchos::RCP<Epetra_MultiVector>       weighted_basis_mean,
                          Teuchos::RCP<Epetra_MultiVector>       kernel_c           ,
                          bool                                   project);

      virtual void Solve();

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y);

      virtual PreconditionerType * Preconditioner() { return NULL; }

      bool IsFactored() { return factored_; }

    private:

      std::string solvertype_;

      //! flag indicating whether matrix was factored before
      bool                                  factored_;

      //! a linear problem wrapper class used by Trilinos and for scaling of the system
      Teuchos::RCP<Epetra_LinearProblem>             lp_;

      //! initial guess and solution
      Teuchos::RCP<Epetra_MultiVector>               x_;

      //! right hand side vector
      Teuchos::RCP<Epetra_MultiVector>               b_;

      //! system of equations
      Teuchos::RCP<Epetra_Operator>  A_;

      //! an abstract amesos solver that can be any of the amesos concrete implementations
      Teuchos::RCP<Amesos_BaseSolver>                amesos_;

      //! reindex linear problem for amesos
      Teuchos::RCP<EpetraExt::LinearProblem_Reindex> reindexer_;
    };
  }
}


#endif /* SOLVER_DIRECTSOLVER_H_ */
