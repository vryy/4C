/*
 * solver_krylovsolver.H
 *
 *  Created on: Jul 4, 2011
 *      Author: wiesner
 */

#ifndef SOLVER_KRYLOVSOLVER_H_
#define SOLVER_KRYLOVSOLVER_H_

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AztecOO.h"
#include "AztecOO_StatusTestResNorm.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestMaxIters.h"

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

          //! an aztex maxiter status test
          Teuchos::RCP<AztecOO_StatusTestMaxIters>       aztest_maxiter_;
          //! an aztec 2-norm convergence test
          Teuchos::RCP<AztecOO_StatusTestResNorm>        aztest_norm2_;
          //! an aztec inf-norm convergence test
          Teuchos::RCP<AztecOO_StatusTestResNorm>        aztest_norminf_;
          // ! combination of the above tests
          Teuchos::RCP<AztecOO_StatusTestCombo>          aztest_combo1_;
          Teuchos::RCP<AztecOO_StatusTestCombo>          aztest_combo2_;

#ifdef WRITEOUTSTATISTICS
          double dtimeprecondsetup_;
#endif
    };
  }
}

#endif /* SOLVER_KRYLOVSOLVER_H_ */
