/*!----------------------------------------------------------------------
\file solver_solvertype.H

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_SOLVERTYPE_H_
#define SOLVER_SOLVERTYPE_H_

#include "../drt_lib/drt_dserror.H"

namespace LINALG
{
  class KrylovProjector;
  namespace SOLVER
  {
    class PreconditionerType;

    /// linear solver type base class
    class SolverType
    {
    public:

      virtual ~SolverType() {}

      /// solve the linear system
      virtual void Setup( Teuchos::RCP<Epetra_Operator>     matrix             ,
                          Teuchos::RCP<Epetra_MultiVector>       x             ,
                          Teuchos::RCP<Epetra_MultiVector>       b             ,
                          bool                             refactor            ,
                          bool                             reset               ,
                          Teuchos::RCP<LINALG::KrylovProjector> projector) = 0;

      virtual int Solve() = 0;

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) = 0;

      /// preconditioner used within this solver
      virtual PreconditionerType * Preconditioner() = 0;

      /// return number of iterations performed by solver
      virtual int getNumIters() const { dserror("Not implemented in base class!"); return -1; };
    };
  }
}


#endif /* SOLVER_SOLVERTYPE_H_ */
