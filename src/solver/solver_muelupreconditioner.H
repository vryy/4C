/*!----------------------------------------------------------------------
\file solver_mlpreconditioner.H

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
Created on: Dec 1, 2011
</pre>
*----------------------------------------------------------------------*/

#ifndef SOLVER_MUELUPRECONDITIONER_H_
#define SOLVER_MUELUPRECONDITIONER_H_

#ifdef HAVE_MueLu

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    /// ml preconditioners
    /*!
      Set of single-matrix algebraic multi-grid preconditioners.
     */
    class MueLuPreconditioner : public PreconditionerType
    {
    public:

      MueLuPreconditioner( FILE * outfile, Teuchos::ParameterList & mllist );

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*P_; }

      virtual void Print( std::ostream & stream ) { stream << "MueLu"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "MueLu Paramters"; }
    private:

      Teuchos::ParameterList & mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_CrsMatrix>                 Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator>                  P_;

    };
  }
}

#endif // HAVE_MueLu

#endif /* SOLVER_MUELUPRECONDITIONER_H_ */
