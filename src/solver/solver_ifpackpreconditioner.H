/*----------------------------------------------------------------------------*/
/*! \file

\brief LINALG::SOLVER wrapper around Trilinos' IFPACK preconditioner

\level 0

*/
/*----------------------------------------------------------------------------*/
#ifndef SOLVER_IFPACKPRECONDITIONER_H_
#define SOLVER_IFPACKPRECONDITIONER_H_

#include "Ifpack.h"

#include "solver_preconditionertype.H"

namespace LINALG::SOLVER
{
  /*! \brief  IFPACK preconditioners
   *
   *  Set of standard single-matrix preconditioners.
   */
  class IFPACKPreconditioner : public SOLVER::PreconditionerType
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
    void Print(std::ostream& stream) override { stream << "IFPACK"; }

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
}  // namespace LINALG::SOLVER

#endif
