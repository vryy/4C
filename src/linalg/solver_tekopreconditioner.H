/*
 * solver_tekopreconditioner.H
 *
 *  Created on: Jul 11, 2011
 *      Author: wiesner
 */

#ifndef SOLVER_TEKOPRECONDITIONER_H_
#define SOLVER_TEKOPRECONDITIONER_H_

#ifdef TRILINOS_DEV

#include "linalg_blocksparsematrix.H"

// Teko specific includes
#include <Teko_EpetraBlockPreconditioner.hpp>

#include "teko_baciepetraoperatorwrapper.H"

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

      TekoPreconditioner( FILE * outfile, Teuchos::ParameterList & tekolist);

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*prec_; }

      virtual void Print( std::ostream & stream ) { stream << "TEKO"; }

      Teuchos::ParameterList& Params() {return tekolist_; };

    private:

      Teuchos::ParameterList & tekolist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<LINALG::BlockSparseMatrixBase>       Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Teko::Epetra::EpetraBlockPreconditioner> prec_;
    };
  }
}

#endif /* TRILINOS_DEV */
#endif /* SOLVER_TEKOPRECONDITIONER_H_ */

