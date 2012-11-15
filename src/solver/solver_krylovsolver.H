/*!----------------------------------------------------------------------
\file solver_krylovsolver.H

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/

#ifndef SOLVER_KRYLOVSOLVER_H_
#define SOLVER_KRYLOVSOLVER_H_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "solver_solvertype.H"
#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {

    /// kryolv subspace linear solvers (aztec) with right-side preconditioning
    class KrylovSolver : public SolverType
    {
      public:

      KrylovSolver( const Epetra_Comm & comm, Teuchos::ParameterList & params, FILE * outfile );

      virtual ~KrylovSolver();

      virtual void Setup( Teuchos::RCP<Epetra_Operator>     matrix             ,
          Teuchos::RCP<Epetra_MultiVector>       x             ,
          Teuchos::RCP<Epetra_MultiVector>       b             ,
          bool                             refactor            ,
          bool                             reset               ,
          Teuchos::RCP<Epetra_MultiVector>  weighted_basis_mean,
          Teuchos::RCP<Epetra_MultiVector>  kernel_c           ,
          bool                             project) = 0;

      virtual void Solve() = 0;

      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y);

      virtual PreconditionerType * Preconditioner() { return &*preconditioner_; }

      int Ncall() { return ncall_; }

      Teuchos::ParameterList & Params() const { return params_; }

      protected:

      //! function for creating preconditioner object
      void CreatePreconditioner( Teuchos::ParameterList & azlist,
          bool isCrsMatrix,
          Teuchos::RCP<Epetra_MultiVector> weighted_basis_mean,
          Teuchos::RCP<Epetra_MultiVector> kernel_c,
          bool project );

          //! a communicator
          const Epetra_Comm&                    comm_;

          //! (internal) parameter list
          Teuchos::ParameterList&                  params_;

          //! file to write output to (proc 0 only, can be NULL on input)
          FILE*                                 outfile_;

          //! initial guess and solution
          Teuchos::RCP<Epetra_MultiVector>                    x_;

          //! right hand side vector
          Teuchos::RCP<Epetra_MultiVector>                    b_;

          //! system of equations
          Teuchos::RCP<Epetra_Operator>  A_;

          //! counting how many times matrix was solved between resets
          int                                   ncall_;

          //! preconditioner object
          Teuchos::RCP<PreconditionerType> preconditioner_;

#ifdef WRITEOUTSTATISTICS
          double dtimeprecondsetup_;
#endif
    };
  }
}

#endif /* SOLVER_KRYLOVSOLVER_H_ */
