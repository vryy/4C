/*!----------------------------------------------------------------------

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
*----------------------------------------------------------------------*/
#ifndef SOLVER_KRYLOVPROJECTIONPRECONDITIONER_H_
#define SOLVER_KRYLOVPROJECTIONPRECONDITIONER_H_

#include "solver_preconditionertype.H"

// forward declarations
class Epetra_Operator;
class Epetra_MultiVector;

namespace LINALG
{
  class LinalgProjectedOperator;
  class KrylovProjector;

  namespace SOLVER
  {
    /// krylov projection for undefined pressure value in incompressible fluids
    /*!
      This is not a preconditioner in a mathematical sense, but it fits the
      software framework nicely.

      A "real" preconditioner is wrapped.
     */
    class KrylovProjectionPreconditioner : public PreconditionerType
    {
     public:
      KrylovProjectionPreconditioner(FILE* outfile, Teuchos::RCP<PreconditionerType> preconditioner,
          Teuchos::RCP<LINALG::KrylovProjector> projector);

      // virtual Epetra_LinearProblem & LinearProblem() { return preconditioner_->LinearProblem(); }

      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      /// return the projecting operator
      virtual Epetra_Operator* PrecOperator() const { return &*P_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print(std::ostream& stream)
      {
        stream << "krylov projection(";
        preconditioner_->Print(stream);
        stream << ")";
      }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "unknown"; }

     private:
      Teuchos::RCP<PreconditionerType> preconditioner_;

      /// Peter's projector object that does the actual work
      Teuchos::RCP<LINALG::KrylovProjector> projector_;

      /// linear operator that calls a "real" preconditioning operator and does
      /// a projection afterwards.
      Teuchos::RCP<LINALG::LinalgProjectedOperator> A_;

      Teuchos::RCP<Epetra_Operator> P_;
    };

  }  // namespace SOLVER
}  // namespace LINALG


#endif /* SOLVER_KRYLOVPROJECTIONPRECONDITIONER_H_ */
