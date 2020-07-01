/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of interface to AztecOO solver package

\level 0

*/
/*---------------------------------------------------------------------*/
#ifndef SOLVER_AZTECSOLVER_H_
#define SOLVER_AZTECSOLVER_H_

#include <MueLu_PermutationFactory_fwd.hpp>
#include <MueLu_FactoryBase_fwd.hpp>
#include <MueLu_Level_fwd.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

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
      //! Constructor
      AztecSolver(const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile);

      //! Destructor
      virtual ~AztecSolver();

      /*! \brief Setup the solver object
       *
       * @param matrix Matrix of the linear system
       * @param x Solution vector of the linear system
       * @param b Right-hand side vector of the linear system
       * @param refactor Boolean flag to enforce a refactorization of the matrix
       * @param reset Boolean flag to enforce a full reset of the solver object
       * @param projector Krylov projector
       */
      virtual void Setup(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
          Teuchos::RCP<Epetra_MultiVector> b, const bool refactor, const bool reset,
          Teuchos::RCP<LINALG::KrylovProjector> projector);

      //! Actual call to the underlying Aztec solver
      virtual int Solve();

      //! return number of iterations
      int getNumIters() const { return numiters_; };

     private:
      //! the Krylov projector
      Teuchos::RCP<LINALG::KrylovProjector> projector_;
      //! an aztec maxiter status test
      Teuchos::RCP<AztecOO_StatusTestMaxIters> aztest_maxiter_;
      //! an aztec projected residual 2-norm convergence test
      Teuchos::RCP<AztecOO_StatusTestProjResNorm> aztest_norm2_;
      //! a test combination (for norm2 and maxiter)
      Teuchos::RCP<AztecOO_StatusTestCombo> aztest_combo_;
      //! number of iterations
      int numiters_;
    };
  }  // namespace SOLVER
}  // namespace LINALG


#endif /* SOLVER_AZTECSOLVER_H_ */
