/*!----------------------------------------------------------------------
\file solver_belossolver.H

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_BELOSSOLVER_H_
#define SOLVER_BELOSSOLVER_H_

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "solver_solvertype.H"
#include "solver_preconditionertype.H"

namespace LINALG
{

  namespace SOLVER
  {
    /// krylov subspace linear solvers (belos) with right-side preconditioning
    class BelosSolver : public SolverType
    {
      public:
        BelosSolver( const Epetra_Comm& comm, Teuchos::ParameterList & params, FILE* outfile);
        virtual ~BelosSolver();

        virtual void Setup( Teuchos::RCP<Epetra_Operator>     matrix             ,
            Teuchos::RCP<Epetra_Vector>       x                  ,
            Teuchos::RCP<Epetra_Vector>       b                  ,
            bool                             refactor           ,
            bool                             reset              ,
            Teuchos::RCP<Epetra_MultiVector>  weighted_basis_mean,
            Teuchos::RCP<Epetra_MultiVector>  kernel_c           ,
            bool                             project);

        virtual void Solve();

        virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y);

        virtual PreconditionerType * Preconditioner() { return &*preconditioner_; }

        int Ncall() { return ncall_; }

        Teuchos::ParameterList & Params() const { return params_; }

      private:

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
        Teuchos::RCP<Epetra_Vector>                    x_;

        //! right hand side vector
        Teuchos::RCP<Epetra_Vector>                    b_;

        //! system of equations
        Teuchos::RCP<Epetra_Operator>  A_;

        //! counting how many times matrix was solved between resets
        int                                   ncall_;

        Teuchos::RCP<PreconditionerType> preconditioner_;

      };
  }
}


#endif /* SOLVER_BELOSSOLVER_H_ */
