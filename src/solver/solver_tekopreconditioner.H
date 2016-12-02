/*!----------------------------------------------------------------------
\file solver_tekopreconditioner.H

\brief Preconditioner class using Trilinos TEKOS framework

\level 2

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de

Created on: Jul 13, 2011
*----------------------------------------------------------------------*/

#ifndef SOLVER_TEKOPRECONDITIONER_H_
#define SOLVER_TEKOPRECONDITIONER_H_

#ifdef HAVE_TEKO

#include "../linalg/linalg_blocksparsematrix.H"

// Teko specific includes
#include <Teko_EpetraBlockPreconditioner.hpp>

#include "teko/teko_baciepetraoperatorwrapper.H"

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    /// Teko preconditioners
    /*!
      Set of standard block-matrix preconditioners

      Teko only needs the parameters for the block inverses in the parameter list
      as sublists "Inverse1", "Inverse2",...
      From the parameter list the Teko lists are automatically constructed, no BACI
      SOLVER objects needed!
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

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return prec_; }

      virtual void Print( std::ostream & stream )
      {
        stream << "TEKO " << type_ << " preconditioner";
      }

      Teuchos::ParameterList& Params() {return tekolist_; };

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "Teko Parameters"; }
    private:

      Teuchos::ParameterList & tekolist_;

      std::string type_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<LINALG::BlockSparseMatrixBase>       Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Teko::Epetra::EpetraBlockPreconditioner> prec_;
    };
  }
}

#endif /* HAVE_TEKO */
#endif /* SOLVER_TEKOPRECONDITIONER_H_ */
