/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
*/
/*----------------------------------------------------------------------*/

#ifndef SOLVER_MLPRECONDITIONER_H_
#define SOLVER_MLPRECONDITIONER_H_

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    /// ml preconditioners
    /*!
      Set of single-matrix algebraic multi-grid preconditioners.
     */
    class MLPreconditioner : public PreconditionerType
    {
     public:
      MLPreconditioner(FILE* outfile, Teuchos::ParameterList& mllist);

      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual Epetra_Operator* PrecOperator() const { return &*P_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print(std::ostream& stream) { stream << "ML"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "ML Parameters"; }

     private:
      Teuchos::ParameterList& mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_RowMatrix> Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator> P_;
    };
  }  // namespace SOLVER
}  // namespace LINALG

#endif /* SOLVER_MLPRECONDITIONER_H_ */
