/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_PRECONDITIONER_BLOCK_H
#define BACI_LINEAR_SOLVER_PRECONDITIONER_BLOCK_H

#include "baci_linear_solver_preconditioner_type.H"

BACI_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
{
  /// SIMPLE(R) block preconditioner
  /*!
    Block preconditioners assume the Epetra_Operator to be a
    CORE::LINALG::BlockSparseMatrix.
   */
  class SimplePreconditioner : public PreconditionerType
  {
   public:
    SimplePreconditioner(Teuchos::ParameterList& params);

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override { return P_; }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    std::string getParameterListName() const override { return "CheapSIMPLE Parameters"; }

   private:
    Teuchos::ParameterList& params_;
    Teuchos::RCP<Epetra_Operator> P_;
  };

  /// General purpose block gauss-seidel preconditioner
  /*!
    2x2 block preconditioner
   */
  class BGSPreconditioner : public PreconditionerType
  {
   public:
    BGSPreconditioner(Teuchos::ParameterList& params, Teuchos::ParameterList& bgslist);

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override { return P_; }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    std::string getParameterListName() const override { return "BGS Parameters"; }

   private:
    Teuchos::ParameterList& params_;
    Teuchos::ParameterList& bgslist_;
    Teuchos::RCP<Epetra_Operator> P_;
  };
}  // namespace CORE::LINEAR_SOLVER

BACI_NAMESPACE_CLOSE

#endif
