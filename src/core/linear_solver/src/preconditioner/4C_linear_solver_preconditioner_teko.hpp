/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for Teko block preconditioner

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_TEKO_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_TEKO_HPP

#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linear_solver_preconditioner_type.hpp"

#include <MueLu_UseDefaultTypes.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinearSolver
{
  /*! \brief Set of standard block-matrix preconditioners

    Teko only needs the parameters for the block inverses in the parameter list
    as sublists "Inverse1", "Inverse2",...
    From the parameter list the Teko lists are automatically constructed, no 4C
    SOLVER objects needed!
   */
  class TekoPreconditioner : public PreconditionerTypeBase
  {
   public:
    TekoPreconditioner(Teuchos::ParameterList& tekolist);

    void setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> prec_operator() const override { return p_; }

   private:
    Teuchos::ParameterList& tekolist_;

    //! system of equations used for preconditioning used by P_ only
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> pmatrix_;

    //! preconditioner
    Teuchos::RCP<Epetra_Operator> p_;
  };
}  // namespace Core::LinearSolver

FOUR_C_NAMESPACE_CLOSE

#endif
