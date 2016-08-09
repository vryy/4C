/*!----------------------------------------------------------------------
\file solver_pointpreconditioner.H

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_POINTPRECONDITIONER_H_
#define SOLVER_POINTPRECONDITIONER_H_

#include "../linalg/linalg_downwindmatrix.H"
#include "../linalg/linalg_none_operator.H"

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    class NonePreconditioner : public PreconditionerType
    {
      public:

      NonePreconditioner( FILE * outfile, Teuchos::ParameterList & list );

      virtual void Setup( bool create,
          Epetra_Operator * matrix,
          Epetra_MultiVector * x,
          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*prec_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return prec_; }

      virtual void Print( std::ostream & stream ) { stream << "None"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "unknown"; }

      private:

      Teuchos::RCP<LINALG::NoneOperator> prec_;
    };



    /// downwinding, a preconditioner-wrapper
    /*!
    Downwinding does a matrix-reordering that is supposed to speed-up the
    solution process. The actual solution is performed using another
    preconditioner. DWindPreconditioner is just a wrapper that presents a
    permuted matrix to the actual preconditioner (and solver).

    \note DWindPreconditioner is a preconditioner in the solver framework, not
    in a mathematical sense.
     */
    class DWindPreconditioner : public PreconditionerType
    {
      public:

        DWindPreconditioner( FILE * outfile,
            Teuchos::RCP<PreconditionerType> preconditioner,
            Teuchos::ParameterList & azlist );

        virtual void Setup( bool create,
            Epetra_Operator * matrix,
            Epetra_MultiVector * x,
            Epetra_MultiVector * b );

        /// undo reordering of lhs, don't care for rhs or matrix
        virtual void Finish( Epetra_Operator * matrix,
            Epetra_MultiVector * x,
            Epetra_MultiVector * b );

        /// return the wrapped preconditioner, no work needs to be done after the
        /// permutation is performed
        virtual Epetra_Operator * PrecOperator() const { return preconditioner_->PrecOperator(); }

        /// linear operator used for preconditioning
        virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return preconditioner_->PrecOperatorRCP(); }

        virtual void Print( std::ostream & stream )
        {
          stream << "DWind(";
          preconditioner_->Print( stream );
          stream << ")";
        }

        /// return name of sublist in paramterlist which contains parameters for preconditioner
        virtual const std::string getParameterListName() const { return "unknown"; }

      private:

        Teuchos::RCP<SOLVER::PreconditionerType> preconditioner_;

        Teuchos::ParameterList & azlist_;

        //! downwinding reindexer for downwind smoothers
        Teuchos::RCP<LINALG::DownwindMatrix>            dwind_;

        Teuchos::RCP<Epetra_CrsMatrix>     dwA_;
        Teuchos::RCP<Epetra_MultiVector>   dwx_;
        Teuchos::RCP<Epetra_MultiVector>   dwb_;
    };

    /// inf-norm scaling fits the preconditioner framework perfecty
    /*!
    Modifies the underlying matrix and needs to unscale the result
    afterwards. Can be combined with any single-matrix preconditioner.
     */
    class InfNormPreconditioner : public PreconditionerType
    {
      public:

      InfNormPreconditioner( Teuchos::RCP<PreconditionerType> preconditioner );

      virtual Epetra_LinearProblem & LinearProblem() { return preconditioner_->LinearProblem(); }

      virtual void Setup( bool create,
          Epetra_Operator * matrix,
          Epetra_MultiVector * x,
          Epetra_MultiVector * b );

      virtual void Finish( Epetra_Operator * matrix,
          Epetra_MultiVector * x,
          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return preconditioner_->PrecOperator(); }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const  { return preconditioner_->PrecOperatorRCP(); }

      virtual void Print( std::ostream & stream )
      {
        stream << "InfNorm(";
        preconditioner_->Print( stream );
        stream << ")";
      }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "unknown"; }

      private:
      Teuchos::RCP<PreconditionerType> preconditioner_;
      Teuchos::RCP<Epetra_Vector> rowsum_;
      Teuchos::RCP<Epetra_Vector> colsum_;
    };

    /// diagonal scaling fits the preconditioner framework perfecty
    /*!
    Modifies the underlying matrix and needs to unscale the result
    afterwards. Can be combined with any single-matrix preconditioner.
     */
    class SymDiagPreconditioner : public PreconditionerType
    {
      public:

      SymDiagPreconditioner( Teuchos::RCP<PreconditionerType> preconditioner );

      virtual Epetra_LinearProblem & LinearProblem() { return preconditioner_->LinearProblem(); }

      virtual void Setup( bool create,
          Epetra_Operator * matrix,
          Epetra_MultiVector * x,
          Epetra_MultiVector * b );

      virtual void Finish( Epetra_Operator * matrix,
          Epetra_MultiVector * x,
          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return preconditioner_->PrecOperator(); }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const  { return preconditioner_->PrecOperatorRCP(); }

      virtual void Print( std::ostream & stream )
      {
        stream << "SymDiag(";
        preconditioner_->Print( stream );
        stream << ")";
      }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "unknown"; }

      private:

      /// embedded preconditioner
      Teuchos::RCP<SOLVER::PreconditionerType> preconditioner_;

      /// matrix diagonal
      Teuchos::RCP<Epetra_Vector> diag_;
    };

  }
}


#endif /* SOLVER_POINTPRECONDITIONER_H_ */
