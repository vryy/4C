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

#include <Amesos_Klu.h>
#ifndef HAVENOT_UMFPACK
#include <Amesos_Umfpack.h>
#endif
#include <Amesos_Lapack.h>

// Trilinos is configured with SuperLUDIST only in the parallel version (not on altix)
#ifdef PARALLEL
#ifndef HAVENOT_SUPERLU
#include <Amesos_Superludist.h>
#endif
#endif

#include <EpetraExt_Reindex_LinearProblem.h>

#include "solver_solvertype.H"

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

      virtual void Setup( RCP<Epetra_Operator>     matrix             ,
                          RCP<Epetra_Vector>       x                  ,
                          RCP<Epetra_Vector>       b                  ,
                          bool                             refactor           ,
                          bool                             reset              ,
                          RCP<Epetra_MultiVector>  weighted_basis_mean,
                          RCP<Epetra_MultiVector>  kernel_c           ,
                          bool                             project);

      virtual void Solve();

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y);

      virtual PreconditionerType * Preconditioner() { return NULL; }

      bool IsFactored() { return factored_; }

    private:

      std::string solvertype_;

      //! flag indicating whether matrix was factored before
      bool                                  factored_;

      //! a linear problem wrapper class used by Trilinos and for scaling of the system
      RCP<Epetra_LinearProblem>             lp_;

      //! initial guess and solution
      RCP<Epetra_Vector>                    x_;

      //! right hand side vector
      RCP<Epetra_Vector>                    b_;

      //! system of equations
      RCP<Epetra_Operator>  A_;

      //! an abstract amesos solver that can be any of the amesos concrete implementations
      RCP<Amesos_BaseSolver>                amesos_;

      //! reindex linear problem for amesos
      RCP<EpetraExt::LinearProblem_Reindex> reindexer_;
    };
  }
}


#endif /* SOLVER_DIRECTSOLVER_H_ */
