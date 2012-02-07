/*!----------------------------------------------------------------------
\file solver_aztecsolver.H

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_AZTECSOLVER_H_
#define SOLVER_AZTECSOLVER_H_

#include "solver_krylovsolver.H"

// Aztec forward declarations
class AztecOO_StatusTestResNorm;
class AztecOO_StatusTestCombo;
class AztecOO_StatusTestMaxIters;

namespace LINALG
{
  namespace SOLVER
  {
    /// kryolv subspace linear solvers (aztec) with right-side preconditioning
    class AztecSolver : public KrylovSolver
    {
      public:

      AztecSolver( const Epetra_Comm & comm, Teuchos::ParameterList & params, FILE * outfile );

      virtual ~AztecSolver();

      virtual void Setup( Teuchos::RCP<Epetra_Operator>     matrix             ,
          Teuchos::RCP<Epetra_MultiVector>       x             ,
          Teuchos::RCP<Epetra_MultiVector>       b             ,
          bool                             refactor            ,
          bool                             reset               ,
          Teuchos::RCP<Epetra_MultiVector>  weighted_basis_mean,
          Teuchos::RCP<Epetra_MultiVector>  kernel_c           ,
          bool                             project);

      virtual void Solve();

      private:

          //! an aztex maxiter status test
          Teuchos::RCP<AztecOO_StatusTestMaxIters>       aztest_maxiter_;
          //! an aztec 2-norm convergence test
          Teuchos::RCP<AztecOO_StatusTestResNorm>        aztest_norm2_;
          //! an aztec inf-norm convergence test
          Teuchos::RCP<AztecOO_StatusTestResNorm>        aztest_norminf_;
          // ! combination of the above tests
          Teuchos::RCP<AztecOO_StatusTestCombo>          aztest_combo1_;
          Teuchos::RCP<AztecOO_StatusTestCombo>          aztest_combo2_;

    };
  }
}


#endif /* SOLVER_AZTECSOLVER_H_ */
