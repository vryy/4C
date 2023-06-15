/*----------------------------------------------------------------------------*/
/*! \file

\brief LINALG::SOLVER wrapper around Trilinos' IFPACK preconditioner

\level 0

*/
/*----------------------------------------------------------------------------*/
#ifndef LINEAR_SOLVER_PRECONDITIONER_IFPACK_H
#define LINEAR_SOLVER_PRECONDITIONER_IFPACK_H

#include <Ifpack.h>

#include "linear_solver_preconditioner_type.H"

namespace CORE::LINEAR_SOLVER
{
  /*! \brief  IFPACK preconditioners
   *
   *  Set of standard single-matrix preconditioners.
   */
  class IFPACKPreconditioner : public LINEAR_SOLVER::PreconditionerType
  {
   public:
    //! Constructor (empty)
    IFPACKPreconditioner(
        FILE* outfile, Teuchos::ParameterList& ifpacklist, Teuchos::ParameterList& azlist);

    //! Setup
    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    //! Access to preconditioner operator
    Epetra_Operator* PrecOperator() const override { return &*prec_; }

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const override { return prec_; }

    //! Print name of preconditioner
    void Print(std::ostream& stream) const override { stream << "IFPACK"; }

    //! return name of sublist in parameter list which contains parameters for preconditioner
    const std::string getParameterListName() const override { return "IFPACK Parameters"; }

   private:
    //! IFPACK parameter list
    Teuchos::ParameterList& ifpacklist_;

    //! solver parameter list
    Teuchos::ParameterList& solverlist_;

    //! system of equations used for preconditioning used by P_ only
    Teuchos::RCP<Epetra_RowMatrix> Pmatrix_;

    //! preconditioner
    Teuchos::RCP<Ifpack_Preconditioner> prec_;

  };  // class IFPACKPreconditioner
}  // namespace CORE::LINEAR_SOLVER

#endif
