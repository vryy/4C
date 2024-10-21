// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_IFPACK_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_IFPACK_HPP

#include "4C_config.hpp"

#include "4C_linear_solver_preconditioner_type.hpp"

#include <Ifpack.h>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinearSolver
{
  /*! \brief  IFPACK preconditioners
   *
   *  Set of standard single-matrix preconditioners.
   */
  class IFPACKPreconditioner : public LinearSolver::PreconditionerTypeBase
  {
   public:
    //! Constructor (empty)
    IFPACKPreconditioner(Teuchos::ParameterList& ifpacklist, Teuchos::ParameterList& solverlist);

    //! Setup
    void setup(bool create, Epetra_Operator* matrix, Core::LinAlg::MultiVector<double>* x,
        Core::LinAlg::MultiVector<double>* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> prec_operator() const override { return prec_; }

   private:
    //! IFPACK parameter list
    Teuchos::ParameterList& ifpacklist_;

    //! solver parameter list
    Teuchos::ParameterList& solverlist_;

    //! system of equations used for preconditioning used by P_ only
    Teuchos::RCP<Epetra_RowMatrix> pmatrix_;

    //! preconditioner
    Teuchos::RCP<Ifpack_Preconditioner> prec_;

  };  // class IFPACKPreconditioner
}  // namespace Core::LinearSolver

FOUR_C_NAMESPACE_CLOSE

#endif
