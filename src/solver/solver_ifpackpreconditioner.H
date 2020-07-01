/*----------------------------------------------------------------------------*/
/*! \file

\level 0

\brief LINALG::SOLVER wrapper around Trilinos' IFPACK preconditioner
*/
/*----------------------------------------------------------------------------*/
#ifndef SOLVER_IFPACKPRECONDITIONER_H_
#define SOLVER_IFPACKPRECONDITIONER_H_

#include "Ifpack.h"

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    /*! \brief  IFPACK preconditioners
     *
     *  Set of standard single-matrix preconditioners.
     */
    class IFPACKPreconditioner : public SOLVER::PreconditionerType
    {
     public:
      //! Constructor (empty)
      IFPACKPreconditioner(FILE* outfile,
          Teuchos::ParameterList& ifpacklist,  ///< parameter list for IFPACK
          Teuchos::ParameterList& azlist       ///< parameter list for Aztec
      );

      //! Setup
      virtual void Setup(bool create,  ///< enforce rebuiling of the preconditioner
          Epetra_Operator* matrix,     ///< matrix
          Epetra_MultiVector* x,       ///< solution vector
          Epetra_MultiVector* b        ///< right-hand side vector
      );

      //! Access to preconditioner operator
      virtual Epetra_Operator* PrecOperator() const { return &*prec_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return prec_; }

      //! Print name of preconditioner
      virtual void Print(std::ostream& stream) { stream << "IFPACK"; }

      //! return name of sublist in parameter list which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "IFPACK Parameters"; }

     private:
      //! IFPACK parameter list
      Teuchos::ParameterList& ifpacklist_;

      //! AZTEC parameter list
      Teuchos::ParameterList& azlist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_RowMatrix> Pmatrix_;

      //! preconditioner
      Teuchos::RCP<Ifpack_Preconditioner> prec_;

    };  // class IFPACKPreconditioner
  }     // namespace SOLVER
}  // namespace LINALG

#endif /* SOLVER_IFPACKPRECONDITIONER_H_ */
