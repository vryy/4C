/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of solver type base class

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_METHOD_HPP
#define FOUR_C_LINEAR_SOLVER_METHOD_HPP

#include "4C_config.hpp"

#include "4C_linalg_krylov_projector.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
{
  /* A very good tutorial and explanation on how to choose your linear solver and the related
   * parameters can be found here:
   * https://de.mathworks.com/help/matlab/math/iterative-methods-for-linear-systems.html
   *
   * The available methods in 4C are very similar to the ones available in MATLAB:
   *
   * For small problems e.g. smaller than 50.000 global degrees of freedom, choose a direct solver
   * (UMFPACK is very popular and also used as direct solver in MATLAB, you just don't see it).
   *
   * For bigger problems use an iterative method in combination with a preconditioner. A popular
   * choice for symmetric systems is ICC + CG, for non-symmetric problems ILU + GMRES. Also try to
   * use Belos as your iterative solver package of choice!
   *
   * For really big problems use a multigrid preconditioner in combination with an iterative
   * solver. This ensures scalability and realistic computation times. Also try to use MueLU as
   * your multigrid package of choice!
   *
   * In 4C equilibration and reordering is also available, but not for everything yet. An
   * experimental approach on permuting the system matrix to obtain better conditioning and
   * faster solution times is given by the permuation strategy.
   *
   */

  //! Available solvers in the Amesos package and Iterative methods
  enum class SolverType
  {
    umfpack,   ///< Amesos direct solver using UMFPACK
    superlu,   ///< Amesos direct solver using SuperLU_Dist
    belos,     ///< Belos iterative solver
    undefined  ///< undefined solver
  };

  //! Different solvers within the Belos package
  enum class IterativeSolverType
  {
    cg,       ///< cg-solver for symmetric problems
    gmres,    ///< gmres-solver for non-symmetric problems
    bicgstab  ///< bicgstab-solver for non-symmetric problems with small storage
  };

  //! Different preconditioners within the ML, MueLu and Ifpack package
  enum class PreconditionerType
  {
    ilu,           ///< incomplete LU factorization with fill in levels (Ifpack package)
    icc,           ///< incomplete Cholesky factorization for symmetric problems (Ifpack package)
    multigrid_ml,  ///< standard multigrid for structures (ML package, outdated)
    multigrid_ml_fluid,   ///< multigrid for fluid problems (ML package, outdated)
    multigrid_ml_fluid2,  ///< energy optimal multigrid for unsymmetric fluid problems (ML package,
                          ///< outdated)
    multigrid_muelu,      ///< multigrid preconditioner (MueLu package, recommended!)
    multigrid_muelu_fluid,  ///< multigrid preconditioner for blocked fluid problems (MueLu package)
    multigrid_muelu_tsi,    ///< multigrid preconditioner for blocked tsi problems (MueLu package)
    multigrid_muelu_contactsp,  ///< multigrid preconditioner for blocked contact problems in saddle
                                ///< point formulation (MueLu package)
    multigrid_muelu_beamsolid,  ///< multigrid preconditioner for blocked beam solid interaction
                                ///< problems (MueLu package)
    multigrid_muelu_fsi,  ///< multigrid preconditioner for blocked fluid structure interaction
                          ///< problems (MueLu package)
    multigrid_nxn,  ///< multigrid preconditioner for a nxn block matrix (indirectly MueLu package)
    block_gauss_seidel_2x2,  ///< block Gauss-Seidel for 2x2 system (inhouse implementation)
    cheap_simple             ///< CheapSIMPLE for 2x2 systems? (inhouse implementation)
  };

  //! scaling strategies for linear solvers
  enum class ScalingStrategy
  {
    none,       ///< no scaling of the linear system
    symmetric,  ///< symmetric scaling of the linear system
    infnorm     ///< infinity-norm scaling of the linear system
  };

  /// linear solver type base class
  template <class MatrixType, class VectorType>
  class SolverTypeBase
  {
   public:
    //! Destructor
    virtual ~SolverTypeBase() = default;

    /*! \brief Setup the solver object
     *
     * @param A Matrix of the linear system
     * @param x Solution vector of the linear system
     * @param b Right-hand side vector of the linear system
     * @param refactor Boolean flag to enforce a refactorization of the matrix
     * @param reset Boolean flag to enforce a full reset of the solver object
     * @param projector Krylov projector
     */
    virtual void Setup(Teuchos::RCP<MatrixType> A, Teuchos::RCP<VectorType> x,
        Teuchos::RCP<VectorType> b, const bool refactor, const bool reset,
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector) = 0;

    virtual int Solve() = 0;

    /// return number of iterations performed by solver
    virtual int getNumIters() const
    {
      FOUR_C_THROW("Not implemented in base class!");
      return -1;
    };
  };
}  // namespace CORE::LINEAR_SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
