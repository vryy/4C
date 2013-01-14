/*!----------------------------------------------------------------------
\file solver_permutedaztecsolver.H

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
Created on: Dec 7, 2012
</pre>
*----------------------------------------------------------------------*/

#ifndef SOLVER_PERMUTEDAZTECSOLVER_H_
#define SOLVER_PERMUTEDAZTECSOLVER_H_

#ifdef HAVE_MueLu
#ifdef HAVE_EXPERIMENTAL_MueLu

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>

#include <MueLu_PermutationFactory_fwd.hpp>
#include <MueLu_FactoryBase_fwd.hpp>
#include <MueLu_Level_fwd.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int
#include <MueLu_UseShortNames.hpp>

#include "solver_krylovsolver.H"

// Aztec forward declarations
class AztecOO_StatusTestResNorm;
class AztecOO_StatusTestCombo;
class AztecOO_StatusTestMaxIters;

class Epetra_MultiVector;

namespace LINALG
{
  namespace SOLVER
  {
    /// kryolv subspace linear solvers (aztec) with right-side preconditioning
    class PermutedAztecSolver : public KrylovSolver
    {
      public:

      PermutedAztecSolver( const Epetra_Comm & comm, Teuchos::ParameterList & params, FILE * outfile );

      virtual ~PermutedAztecSolver();

      virtual void Setup( Teuchos::RCP<Epetra_Operator>     matrix             ,
          Teuchos::RCP<Epetra_MultiVector>       x             ,
          Teuchos::RCP<Epetra_MultiVector>       b             ,
          bool                             refactor            ,
          bool                             reset               ,
          Teuchos::RCP<Epetra_MultiVector>  weighted_basis_mean,
          Teuchos::RCP<Epetra_MultiVector>  kernel_c           ,
          bool                             project);

      virtual void Solve();

      private:
          void BuildPermutationOperator(const Teuchos::RCP<Epetra_CrsMatrix>& A);
          void PermuteLinearSystem(const Teuchos::RCP<Epetra_CrsMatrix>& A,const Teuchos::RCP<Epetra_MultiVector>& b);
          void ReTransformSolution();

          // parameter xOp only needed for RowMap and DomainMap? (use original matrix A without permutation)
          void PermuteNullSpace(const Teuchos::RCP<Epetra_CrsMatrix>& A);

          Teuchos::RCP<const Epetra_CrsMatrix> GetOperator(const std::string name, const Teuchos::RCP<FactoryBase> & fact);
          Teuchos::RCP<Epetra_CrsMatrix>       GetOperatorNonConst(const std::string name, const Teuchos::RCP<FactoryBase> & fact);

          int CountZerosOnDiagonalEpetra(const Teuchos::RCP<Epetra_CrsMatrix>& A);
          int CountZerosOnDiagonal(const Teuchos::RCP<const Matrix>& xOp);

          void test();
      private:
          bool bPermuteLinearSystem_;

          Teuchos::RCP<MueLu::Level> data_; // data container

          Teuchos::RCP<PermutationFactory> PermFact_; // permutation factory

          //! an aztex maxiter status test
          Teuchos::RCP<AztecOO_StatusTestMaxIters>       aztest_maxiter_;
          //! an aztec 2-norm convergence test
          Teuchos::RCP<AztecOO_StatusTestResNorm>        aztest_norm2_;
          //! an aztec inf-norm convergence test
          Teuchos::RCP<AztecOO_StatusTestResNorm>        aztest_norminf_;
          // ! combination of the above tests
          Teuchos::RCP<AztecOO_StatusTestCombo>          aztest_combo1_;
          Teuchos::RCP<AztecOO_StatusTestCombo>          aztest_combo2_;

    };
  }
}

#endif // HAVE_EXPERIMENTAL_MueLu
#endif // HAVE_MueLu

#endif /* SOLVER_PERMUTEDAZTECSOLVER_H_ */
