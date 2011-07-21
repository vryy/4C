/*
 * solver_tekopreconditioner.H
 *
 *  Created on: Jul 11, 2011
 *      Author: wiesner
 */

#ifndef SOLVER_TEKOPRECONDITIONER_H_
#define SOLVER_TEKOPRECONDITIONER_H_

#ifdef TRILINOS_DEV
#ifdef D_TEKO

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    /// ifpack preconditioners
    /*!
      Set of standard single-matrix preconditioners.
     */
    class TekoPreconditioner : public SOLVER::PreconditionerType
    {
    public:

      TekoPreconditioner( FILE * outfile, Teuchos::ParameterList & tekolist, Teuchos::ParameterList & azlist );

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*prec_; }

      virtual void Print( std::ostream & stream ) { stream << "IFPACK"; }

    private:

      Teuchos::ParameterList & tekolist_;
      Teuchos::ParameterList & azlist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_RowMatrix>                 Pmatrix_;

      /// preconditioner
      //Teuchos::RCP<Ifpack_Preconditioner> prec_;
      Teuchos::RCP<Epetra_Operator> prec_;
    };
  }
}

#endif /* TEKO */
#endif /* TRILINOS_DEV */
#endif /* SOLVER_TEKOPRECONDITIONER_H_ */

