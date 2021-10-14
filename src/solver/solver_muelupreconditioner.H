/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for MueLu preconditioner

\level 1
Created on: Dec 1, 2011
*/
/*----------------------------------------------------------------------*/

#ifndef SOLVER_MUELUPRECONDITIONER_H_
#define SOLVER_MUELUPRECONDITIONER_H_

#include "solver_preconditionertype.H"

// MueLu
#include <MueLu_Hierarchy.hpp>

// Header files for default types must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>

namespace LINALG
{
  namespace SOLVER
  {
    /*! \brief MueLu preconditioner for non-blocked linear systems of equations

     Set of single-matrix algebraic multi-grid preconditioners for non-blocked problems
     */
    class MueLuPreconditioner : public PreconditionerType
    {
     public:
      MueLuPreconditioner(FILE* outfile, Teuchos::ParameterList& mllist);

      /*! \brief Create and compute the preconditioner
       *
       * The saddle-point preconditioner only works for matrices of the type
       * Epetra_CrsMatrix. We check whether the input \c matrix is of proper type
       * and throw an error if not.
       *
       * This routine either re-create the entire preconditioner from scratch or
       * it re-uses the existing preconditioner and only updates the fine level matrix
       * for the Krylov solver.
       *
       * It maintains backward compability to the ML interface!
       *
       * @param create Boolean flag to enforce (re-)creation of the preconditioner
       * @param matrix Epetra_CrsMatrix to be used as input for the preconditioner
       * @param x Solution of the linear system
       * @param b Right-hand side of the linear system
       */
      void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) final;

      //! Get the actual operator used for preconditioning
      Epetra_Operator* PrecOperator() const final { return &*P_; }

      //! linear operator used for preconditioning
      Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const final { return P_; }

      void Print(std::ostream& stream) final { stream << "MueLu"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      const std::string getParameterListName() const final { return "MueLu Parameters"; }

     private:
      //! ML parameter list
      Teuchos::ParameterList& mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_CrsMatrix> Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator> P_;

      /// MueLu hierarchy
      Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H_;
    };  // namespace MueLuPreconditioner
  }     // namespace SOLVER
}  // namespace LINALG

#endif /* SOLVER_MUELUPRECONDITIONER_H_ */
