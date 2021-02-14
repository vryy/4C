/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 2

*/
/*-----------------------------------------------------------------------------*/

#ifndef SOLVER_MUELUCONTACTSPPRECONDITIONER_H_
#define SOLVER_MUELUCONTACTSPPRECONDITIONER_H_

#include "solver_preconditionertype.H"

// MueLu
#include <MueLu_Hierarchy.hpp>

// Header files for default types must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>

// BACI includes
#include "../linalg/linalg_blocksparsematrix.H"

namespace LINALG
{
  namespace SOLVER
  {
    /*! \brief MueLu preconditioner for contact problems in saddle-point formulation

     Set of single-matrix algebraic multi-grid preconditioners for contact/meshtying problems
     */
    class MueLuContactSpPreconditioner : public PreconditionerType
    {
     public:
      MueLuContactSpPreconditioner(FILE* outfile, Teuchos::ParameterList& mllist);

      /*! \brief Create and compute the preconditioner
       *
       * The saddle-point preconditioner only works for matrices of the type
       * BlockSparseMatrix. We assume the structure and block indices to be
       * \f[
       * A = \left[\begin{array}{cc}
       *       A_{11} & A_{12}\\
       *       A_{21} & A_{22}
       *     \end{array}\right]
       * \f]
       * We check whether the input \c matrix is of proper type and throw an
       * error if not.
       *
       * This routine either re-create the entire preconditioner from scratch or
       * it re-uses the existing preconditioner and only updates the fine level matrix
       * for the Krylov solver.
       *
       * @param create Boolean flag to enforce (re-)creation of the preconditioner
       * @param matrix BlockSparseMatrix to be used as input for the preconditioner
       * @param x Solution of the linear system
       * @param b Right-hand side of the linear system
       */
      void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) final;

      //! Get the actual operator used for preconditioning
      Epetra_Operator* PrecOperator() const final { return &*P_; }

      //! Get the actual operator used for preconditioning
      Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const final { return P_; }

      void Print(std::ostream& stream) final { stream << "MueLu"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      const std::string getParameterListName() const final { return "MueLu (Contact) Parameters"; }

     private:
      //! ML parameter list
      Teuchos::ParameterList& mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<BlockSparseMatrixBase> Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator> P_;

      /// MueLu hierarchy
      Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H_;
    };  // class MueLuContactSpPreconditioner
  }     // namespace SOLVER
}  // namespace LINALG

#endif /* SOLVER_MUELUCONTACTSPPRECONDITIONER_H_ */
