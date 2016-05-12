/*!----------------------------------------------------------------------
\file solver_belossolver.H

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_BELOSSOLVER_H_
#define SOLVER_BELOSSOLVER_H_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "solver_krylovsolver.H"
#include "solver_preconditionertype.H"

namespace LINALG
{
  class KrylovProjector;

  namespace SOLVER
  {
    /// krylov subspace linear solvers (belos) with right-side preconditioning
    class BelosSolver : public KrylovSolver
    {
      public:
        BelosSolver(
          const Epetra_Comm& comm,
          Teuchos::ParameterList & params,
          FILE* outfile
          );

        virtual ~BelosSolver();

        virtual void Setup(
          Teuchos::RCP<Epetra_Operator>     matrix             ,
          Teuchos::RCP<Epetra_MultiVector>       x             ,
          Teuchos::RCP<Epetra_MultiVector>       b             ,
          bool                             refactor            ,
          bool                             reset               ,
          Teuchos::RCP<LINALG::KrylovProjector> projector
          );
;
        virtual int Solve();
      };
  }
}


#endif /* SOLVER_BELOSSOLVER_H_ */
