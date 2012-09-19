/*!----------------------------------------------------------------------
\file solver_krylovprojectionpreconditioner.H

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_KRYLOVPROJECTIONPRECONDITIONER_H_
#define SOLVER_KRYLOVPROJECTIONPRECONDITIONER_H_

#include "solver_preconditionertype.H"

class Epetra_Operator;
class Epetra_MultiVector;

namespace LINALG
{
  class LinalgProjectedOperator;
  class KrylovProjector;

  namespace SOLVER
  {
    /// krylov projection for undefined pressure value in incompressible fluids
    /*!
      This is not a preconditioner in a mathematical sense, but it fits the
      software framework nicely.

      A "real" preconditioner is wrapped.
     */
    class KrylovProjectionPreconditioner : public PreconditionerType
    {
    public:

      KrylovProjectionPreconditioner( FILE * outfile,
                                      Teuchos::RCP<PreconditionerType>  preconditioner,
                                      Teuchos::RCP<Epetra_MultiVector>  weighted_basis_mean,
                                      Teuchos::RCP<Epetra_MultiVector>  kernel_c );

      //virtual Epetra_LinearProblem & LinearProblem() { return preconditioner_->LinearProblem(); }

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual void Finish( Epetra_Operator * matrix,
                           Epetra_MultiVector * x,
                           Epetra_MultiVector * b );

      /// return the projecting operator
      virtual Epetra_Operator * PrecOperator() const { return &*P_; }

      virtual void Print( std::ostream & stream )
      {
        stream << "krylov projection(";
        preconditioner_->Print( stream );
        stream << ")";
      }

    private:
      Teuchos::RCP<PreconditionerType> preconditioner_;

      /// Peter's projector object that does the actual work
      Teuchos::RCP<LINALG::KrylovProjector> projector_;
      Teuchos::RCP<Epetra_MultiVector> weighted_basis_mean_;
      Teuchos::RCP<Epetra_MultiVector> kernel_c_;

      /// linear operator that calls a "real" preconditioning operator and does
      /// a projection afterwards.
      Teuchos::RCP<LINALG::LinalgProjectedOperator> A_;

      Teuchos::RCP<Epetra_Operator> P_;
    };

  }
}


#endif /* SOLVER_KRYLOVPROJECTIONPRECONDITIONER_H_ */
