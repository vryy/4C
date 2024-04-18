/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of Baci's interface to Krylov solvers

\level 0

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_METHOD_ITERATIVE_HPP
#define FOUR_C_LINEAR_SOLVER_METHOD_ITERATIVE_HPP

#include "baci_config.hpp"

#include "baci_linear_solver_method.hpp"
#include "baci_linear_solver_preconditioner_type.hpp"
#include "baci_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
{
  //! krylov subspace linear solvers with right-side preconditioning
  template <class MatrixType, class VectorType>
  class IterativeSolver : public SolverType<MatrixType, VectorType>
  {
   public:
    //! Constructor
    IterativeSolver(const Epetra_Comm& comm, Teuchos::ParameterList& params);

    /*! \brief Setup the solver object
     *
     * @param A Matrix of the linear system
     * @param x Solution vector of the linear system
     * @param b Right-hand side vector of the linear system
     * @param refactor Boolean flag to enforce a refactorization of the matrix
     * @param reset Boolean flag to enforce a full reset of the solver object
     * @param projector Krylov projector
     */
    void Setup(Teuchos::RCP<MatrixType> matrix, Teuchos::RCP<VectorType> x,
        Teuchos::RCP<VectorType> b, const bool refactor, const bool reset,
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector) override;

    //! Actual call to the underlying Belos solver
    int Solve() override;

    int Ncall() { return ncall_; }

    //! return number of iterations
    int getNumIters() const override { return numiters_; };

    Teuchos::ParameterList& Params() const { return params_; }

   private:
    //! Get underlying preconditioner object of type PreconditionerType
    const CORE::LINEAR_SOLVER::PreconditionerType& Preconditioner()
    {
      if (preconditioner_ == Teuchos::null)
        FOUR_C_THROW("Preconditioner has to be given for iterative methods!");
      return *preconditioner_;
    }

    /*! \brief Check whether preconditioner will be reused
     *
     * The user can control reuse/recomputation of the preconditioner by setting appropriate input
     * arguments \c reuse and \c reset. In addition, contact mechanics problems perform some
     * additional checks since they require to rebuild the preconditioner when the active set has
     * changed.
     *
     * @param[in] reuse Parameter AZREUSE from parameter list
     * @param reset Force preconditioner to be rebuilt
     * @return Boolean flag to indicate whether preconditioner is reused (\c true) or has to be
     * re-computed (\c false)
     *
     * \sa CheckReuseStatusOfActiveSet()
     */
    bool AllowReusePreconditioner(const int reuse, const bool reset);

    /*! \brief Function for creating preconditioner object
     *
     * @param solverlist liner solver parameter list
     * @param isCrsMatrix Boolean flag to indicate Epetra_CrsMatrix (true) or block matrix (false)
     * @param projector Krylov projector
     */
    Teuchos::RCP<CORE::LINEAR_SOLVER::PreconditionerType> CreatePreconditioner(
        Teuchos::ParameterList& solverlist, const bool isCrsMatrix,
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector);

    //! a communicator
    const Epetra_Comm& comm_;

    //! (internal) parameter lists
    Teuchos::ParameterList& params_;

    //! initial guess and solution
    Teuchos::RCP<VectorType> x_;

    //! right hand side vector
    Teuchos::RCP<VectorType> b_;

    //! system of equations
    Teuchos::RCP<MatrixType> a_;

    //! counting how many times matrix was solved between resets
    int ncall_{0};

    //! number of iterations
    int numiters_{-1};

    //! preconditioner object
    Teuchos::RCP<CORE::LINEAR_SOLVER::PreconditionerType> preconditioner_;

    /*! \brief Check if active set has changed. If yes, enforce to rebuild the preconditioner.
     *
     * We only can reuse the preconditioner if the active set in contact problems has not changed.
     * Therefore, we compare the current map of active DOFs with #activeDofMap_, i.e. the one
     * stored from the previous application of the preconditioner.
     *
     * The comparison is done in multiple stages:
     * 1. Check if number of active DOFs has changed:
     *   - If yes, we need to rebuild the preconditioner.
     *   - If not, then check 2.
     * 2. Compare current and stored map of active DOFs using a PointSameAs() comparison.
     *   - If map has changed, we need to rebuild the preconditioner.
     *
     * \param[in/out] bAllowReuse Boolean flag to indicate reuse (true) or rebuild
     * (false) of the preconditioner
     * \param[in] linSysParams Parameter list with some linear system information
     */
    bool CheckReuseStatusOfActiveSet(const Teuchos::ParameterList& linSysParams);

    /*! \brief Map of active DOFs in structural contact simulations.
     *
     * This is used to check whether preconditioner can be reused or not.
     * We cannot reuse the preconditioner if the active set has changed.
     *
     * \note Member variable only used for structural contact problems with the
     * "contact activeDofMap" parameter in the "Linear System properties"
     * parameter list set.
     *
     * \sa CheckReuseStatusOfActiveSet()
     */
    Teuchos::RCP<Epetra_Map> active_dof_map_;

    //!@}
  };
}  // namespace CORE::LINEAR_SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
