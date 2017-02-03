/*!----------------------------------------------------------------------
\file solver_aztecsolver.H

<pre>
\brief Declaration
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_AZTECSOLVER_H_
#define SOLVER_AZTECSOLVER_H_

#ifdef HAVE_MueLu

#include <MueLu_PermutationFactory_fwd.hpp>
#include <MueLu_FactoryBase_fwd.hpp>
#include <MueLu_Level_fwd.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#endif // HAVE_MueLu

#include "solver_krylovsolver.H"

// Aztec forward declarations
class AztecOO_StatusTestResNorm;
class AztecOO_StatusTestProjResNorm;
class AztecOO_StatusTestCombo;
class AztecOO_StatusTestMaxIters;

namespace LINALG
{
  class KrylovProjector;

  namespace SOLVER
  {
    /// kryolv subspace linear solvers (aztec) with right-side preconditioning
    class AztecSolver : public KrylovSolver
    {
      public:

      AztecSolver( const Epetra_Comm & comm, Teuchos::ParameterList & params, FILE * outfile );

      virtual ~AztecSolver();

      virtual void Setup(
        Teuchos::RCP<Epetra_Operator>     matrix             ,
        Teuchos::RCP<Epetra_MultiVector>       x             ,
        Teuchos::RCP<Epetra_MultiVector>       b             ,
        bool                             refactor            ,
        bool                             reset               ,
        Teuchos::RCP<LINALG::KrylovProjector> projector
        );

      virtual int Solve();

      //! return number of iterations
      int getNumIters() const { return numiters_; };

      private:

      //! the Krylov projector
      Teuchos::RCP<LINALG::KrylovProjector> projector_;
      //! an aztec maxiter status test
      Teuchos::RCP<AztecOO_StatusTestMaxIters>       aztest_maxiter_;
      //! an aztec projected residual 2-norm convergence test
      Teuchos::RCP<AztecOO_StatusTestProjResNorm>    aztest_norm2_;
      //! a test combination (for norm2 and maxiter)
      Teuchos::RCP<AztecOO_StatusTestCombo>          aztest_combo_;
      //! number of iterations
      int numiters_;

    };
  }
}


#endif /* SOLVER_AZTECSOLVER_H_ */
