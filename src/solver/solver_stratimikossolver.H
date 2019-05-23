/*!----------------------------------------------------------------------

\brief Solver class using Trilinos STRATIMIKOS framework

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de

Created on: Jul 13, 2011
*/
/*----------------------------------------------------------------------*/

#ifndef SOLVER_STRATIMIKOSSOLVER_H_
#define SOLVER_STRATIMIKOSSOLVER_H_

#ifdef HAVE_Stratimikos

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Comm.hpp"

#include "solver_solvertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    /// linear solvers (using the Trilinos STRATIMIKOS framework)
    class StratimikosSolver : public SolverType
    {
     public:
      StratimikosSolver(const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile);

      virtual ~StratimikosSolver();

      virtual void Setup(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
          Teuchos::RCP<Epetra_MultiVector> b, bool refactor, bool reset,
          Teuchos::RCP<LINALG::KrylovProjector> projector);

      virtual int Solve();

      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y);

      virtual PreconditionerType* Preconditioner() { return NULL; }

      Teuchos::ParameterList& Params() const { return params_; }

     protected:
      //! a communicator
      const Epetra_Comm& comm_;

      //! (internal) parameter list
      Teuchos::ParameterList& params_;

      //! file to write output to (proc 0 only, can be NULL on input)
      FILE* outfile_;

      //! initial guess and solution
      Teuchos::RCP<Epetra_MultiVector> x_;

      //! right hand side vector
      Teuchos::RCP<Epetra_MultiVector> b_;

      //! system of equations
      Teuchos::RCP<Epetra_Operator> A_;

      //! counting how many times matrix was solved between resets
      int ncall_;
    };
  }  // namespace SOLVER
}  // namespace LINALG

#endif  // #ifdef HAVE_Stratimikos

#endif /* SOLVER_STRATIMIKOSSOLVER_H_ */
