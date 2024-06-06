/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

 \level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_KRYLOVPROJECTION_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_KRYLOVPROJECTION_HPP

#include "4C_config.hpp"

#include "4C_linalg_krylov_projector.hpp"
#include "4C_linalg_projected_operator.hpp"
#include "4C_linear_solver_preconditioner_type.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinearSolver
{
  /// krylov projection for undefined pressure value in incompressible fluids
  /*!
    This is not a preconditioner in a mathematical sense, but it fits the
    software framework nicely. A "real" preconditioner is wrapped.
  */
  class KrylovProjectionPreconditioner : public PreconditionerTypeBase
  {
   public:
    KrylovProjectionPreconditioner(Teuchos::RCP<PreconditionerTypeBase> preconditioner,
        Teuchos::RCP<Core::LinAlg::KrylovProjector> projector);

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override { return p_; }

   private:
    Teuchos::RCP<PreconditionerTypeBase> preconditioner_;

    /// Peter's projector object that does the actual work
    Teuchos::RCP<Core::LinAlg::KrylovProjector> projector_;

    /// linear operator that calls a "real" preconditioning operator and does
    /// a projection afterwards.
    Teuchos::RCP<Core::LinAlg::LinalgProjectedOperator> a_;

    Teuchos::RCP<Epetra_Operator> p_;
  };
}  // namespace Core::LinearSolver

FOUR_C_NAMESPACE_CLOSE

#endif
