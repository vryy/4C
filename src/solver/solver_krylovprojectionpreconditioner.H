/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

 \level 1

*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_KRYLOVPROJECTIONPRECONDITIONER_H
#define SOLVER_KRYLOVPROJECTIONPRECONDITIONER_H

#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>

#include "linalg_krylov_projector.H"
#include "linalg_projected_operator.H"

#include "solver_preconditionertype.H"

namespace LINALG::SOLVER
{
  /// krylov projection for undefined pressure value in incompressible fluids
  /*!
    This is not a preconditioner in a mathematical sense, but it fits the
    software framework nicely. A "real" preconditioner is wrapped.
  */
  class KrylovProjectionPreconditioner : public PreconditionerType
  {
   public:
    KrylovProjectionPreconditioner(FILE* outfile, Teuchos::RCP<PreconditionerType> preconditioner,
        Teuchos::RCP<LINALG::KrylovProjector> projector);

    // virtual Epetra_LinearProblem & LinearProblem() { return preconditioner_->LinearProblem(); }

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) override;

    /// return the projecting operator
    Epetra_Operator* PrecOperator() const override { return &*P_; }

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const override { return P_; }

    void Print(std::ostream& stream) const override
    {
      stream << "krylov projection(";
      preconditioner_->Print(stream);
      stream << ")";
    }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    const std::string getParameterListName() const override { return "unknown"; }

   private:
    Teuchos::RCP<PreconditionerType> preconditioner_;

    /// Peter's projector object that does the actual work
    Teuchos::RCP<LINALG::KrylovProjector> projector_;

    /// linear operator that calls a "real" preconditioning operator and does
    /// a projection afterwards.
    Teuchos::RCP<LINALG::LinalgProjectedOperator> A_;

    Teuchos::RCP<Epetra_Operator> P_;
  };
}  // namespace LINALG::SOLVER

#endif
