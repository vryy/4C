/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

 \level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_KRYLOVPROJECTION_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_KRYLOVPROJECTION_HPP

#include "baci_config.hpp"

#include "baci_linalg_krylov_projector.hpp"
#include "baci_linalg_projected_operator.hpp"
#include "baci_linear_solver_preconditioner_type.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
{
  /// krylov projection for undefined pressure value in incompressible fluids
  /*!
    This is not a preconditioner in a mathematical sense, but it fits the
    software framework nicely. A "real" preconditioner is wrapped.
  */
  class KrylovProjectionPreconditioner : public PreconditionerType
  {
   public:
    KrylovProjectionPreconditioner(Teuchos::RCP<PreconditionerType> preconditioner,
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector);

    // virtual Epetra_LinearProblem & LinearProblem() { return preconditioner_->LinearProblem(); }

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override { return P_; }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    std::string getParameterListName() const override { return "unknown"; }

   private:
    Teuchos::RCP<PreconditionerType> preconditioner_;

    /// Peter's projector object that does the actual work
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector_;

    /// linear operator that calls a "real" preconditioning operator and does
    /// a projection afterwards.
    Teuchos::RCP<CORE::LINALG::LinalgProjectedOperator> A_;

    Teuchos::RCP<Epetra_Operator> P_;
  };
}  // namespace CORE::LINEAR_SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
