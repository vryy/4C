/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_BLOCKPRECONDITIONERS_H_
#define SOLVER_BLOCKPRECONDITIONERS_H_

#include "solver_preconditionertype.H"

namespace LINALG::SOLVER
{
  /// SIMPLE(R) block preconditioner
  /*!
    Block preconditioners assume the Epetra_Operator to be a
    LINALG::BlockSparseMatrix.
   */
  class SimplePreconditioner : public PreconditionerType
  {
   public:
    SimplePreconditioner(FILE* outfile, Teuchos::ParameterList& params);

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    Epetra_Operator* PrecOperator() const override { return &*P_; }

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const override { return P_; }

    void Print(std::ostream& stream) const override { stream << "CheapSimple"; }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    const std::string getParameterListName() const override { return "CheapSIMPLE Parameters"; }

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
    BGSPreconditioner(
        FILE* outfile, Teuchos::ParameterList& params, Teuchos::ParameterList& bgslist);

    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    Epetra_Operator* PrecOperator() const override { return &*P_; }

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const override { return P_; }

    void Print(std::ostream& stream) const override { stream << "BGS"; }

    /// return name of sublist in paramterlist which contains parameters for preconditioner
    const std::string getParameterListName() const override { return "BGS Parameters"; }

   private:
    Teuchos::ParameterList& params_;
    Teuchos::ParameterList& bgslist_;
    Teuchos::RCP<Epetra_Operator> P_;
  };
}  // namespace LINALG::SOLVER

#endif
