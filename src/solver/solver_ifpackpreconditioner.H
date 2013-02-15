/*!----------------------------------------------------------------------
\file solver_ifpackpreconditioner.H

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_IFPACKPRECONDITIONER_H_
#define SOLVER_IFPACKPRECONDITIONER_H_

#include "Ifpack.h"
//#include "Ifpack_AdditiveSchwarz.h"

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    /// ifpack preconditioners
    /*!
      Set of standard single-matrix preconditioners.
     */
    class IFPACKPreconditioner : public SOLVER::PreconditionerType
    {
    public:

      IFPACKPreconditioner( FILE * outfile, Teuchos::ParameterList & ifpacklist, Teuchos::ParameterList & azlist );

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*prec_; }

      virtual void Print( std::ostream & stream ) { stream << "IFPACK"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "IFPACK Parameters"; }
    private:

      Teuchos::ParameterList & ifpacklist_;
      Teuchos::ParameterList & azlist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_RowMatrix>                 Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Ifpack_Preconditioner> prec_;
    };

  }
}

#endif /* SOLVER_IFPACKPRECONDITIONER_H_ */
