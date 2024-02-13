/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of basic solver interface

\level 0

*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_METHOD_LINALG_HPP
#define BACI_LINEAR_SOLVER_METHOD_LINALG_HPP

#include "baci_config.hpp"

#include "baci_inpar_solver.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_LinearProblem;

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class SparseMatrix;
  class SparseOperator;
  class DownwindMatrix;
  class KrylovProjector;
}  // namespace CORE::LINALG

namespace CORE::LINEAR_SOLVER
{
  template <class MatrixType, class VectorType>
  class SolverType;
}

namespace CORE::LINALG
{
  //! parameters to pass to the Solve() call
  struct SolverParams
  {
    //! system should be refactorized
    bool refactor = false;

    //! data from previous solves should be recalculated including preconditioners
    bool reset = false;

    //! Krylov space projector
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector = Teuchos::null;

    //! for adaptivity of the tolerance: tolerance of the nonlinear solver
    double nonlin_tolerance = {};

    //! for adaptivity of the tolerance: current residual of the nonlinear solver
    double nonlin_residual = {};

    //! for adaptivity of the tolerance: factor by which the tolerance of the linear solver should
    //! be stricter than the current residual of the nonlinear solver
    double lin_tol_better = -1.0;

    //! tolerance of linear solver
    double tolerance = -1.0;
  };

  /*!
  \brief A general solver interface to Trilinos solvers and spooles

  - The input of parameters needs to be completely reworked (hiwi job)

  - This class should implement Epetra_Operator in the future

  */
  class Solver
  {
   public:
    /*!
    \brief Destructor

    */
    virtual ~Solver();

    /*!
    \brief Constructor taking a validated input parameter list for Solver

    Creates a solver using the parameters provided by inparams. They are translated
    by #TranslateSolverParameters to the format required by Belos, if @p translate_params is true.
    Otherwise, they need to be provided such that Belos understands them.

    \param inparams (in): input parameter list as provided by GLOBAL::Problem,
                          e.g. GLOBAL::Problem::SolverParams(num)
    \param comm     (in): a reference to a Epetra communicator object

    \param translate_params_to_belos  (in): translate parameters to Belos

    */
    Solver(const Teuchos::ParameterList& inparams, const Epetra_Comm& comm,
        bool translate_params_to_belos = true);

    /*!
    \brief Set-up of stuff common to all constructors

    */
    void Setup();

    //! @name Solve and ApplyInverse methods

    /*!
    \brief Setup system of equations

    \param matrix    (in/out): system of equations
    \param x         (in/out): initial guess on input, solution on output
    \param b         (in/out): right hand side vector
                               if project!=false it might be modified if not orthognal
                               to matrix kernel.
    \param params  (in)    : parameters for the solver. See documentation of SolverParams
    */
    void Setup(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
        Teuchos::RCP<Epetra_MultiVector> b, const SolverParams& params);

    /// solve linear system after setup has been done
    int Solve();

    /*!
    \brief Solve system of equations in one go

    This includes setup. Reuse of preconditioners and factorized systems is
    provided.

    \param matrix    (in/out): system of equations
    \param x         (in/out): initial guess on input, solution on output
    \param b         (in/out): right hand side vector
                               if project!=false it might be modified if not orthognal
                               to matrix kernel.
    \param params  (in)    : parameters for the solver. See documentation of SolverParams
    */
    int Solve(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
        Teuchos::RCP<Epetra_MultiVector> b, const SolverParams& params);

    /*!
    \brief Solve system of equations in one go (NOX)

    Calls the related Solve function. The only difference is the input
    as Epetra_LinearSystem.

    \param linProblem (in/out): linear Problem wrapper. Includes
                                --> system of equations
                                --> initial guess on input, solution on output
                                --> right hand side vector
                                if project!=false it might be modified if not orthognal
                                to matrix kernel.
    \param params  (in)    : parameters for the solver. See documentation of SolverParams
    */
    int NoxSolve(Epetra_LinearProblem& linProblem, const SolverParams& params);

    /*!
    \brief Reset the solver and clear data

    All data is destroyed except the parameter list
    */
    void Reset();

    //! get tolerance from Belos solver
    double GetTolerance() const
    {
      return Params().sublist("Belos Parameters").get<double>("Convergence Tolerance", 1.0e-8);
    }

    /*!
    \brief Adapt tolerance of iterative solver

    Reset the tolerance read from input file. Can only be used after a call to
    AdaptTolerance.

    \note This method works with iterative solvers only - it does nothing for all other
          solvers configured.

    \sa AdaptTolerance
    */
    void ResetTolerance();

    //@}
    //! @name Input of parameters

    /*!
    \brief Translate solver input parameters from input parameter list to
           internal solver parameters list style

    \param inparams (in): input parameter list as provided by GLOBAL::Problem,
                          e.g. GLOBAL::Problem::SolverParams(num) in case of solver for
                          structures and num according to STRUCTURAL DYNAMIC
    \return             : internal parameter list ready to be associated
                          with #params_
    \date 02/07,11/08
    */
    static Teuchos::ParameterList TranslateSolverParameters(const Teuchos::ParameterList& inparams);

    /*!
    \brief Translate BACI dat file parameters

    */
    static Teuchos::ParameterList TranslateBACIToML(
        const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist = nullptr);

    static Teuchos::ParameterList TranslateBACIToMuelu(
        const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist = nullptr);

    static Teuchos::ParameterList TranslateBACIToIfpack(const Teuchos::ParameterList& inparams);

    static Teuchos::ParameterList TranslateBACIToBelos(const Teuchos::ParameterList& inparams);

    /*!
    \brief Add a validated input parameter list as sublist to internal
           parameters

    \param name     (in): name for sublist for #params_ to attach data to
    \param inparams (in): input parameter list as provided by GLOBAL::Problem,
                          e.g. GLOBAL::Problem::FluidPressureSolverParams in case
                          of additional solver for pressure preconditioner
    \author bborn
    \date 11/08
    */
    void PutSolverParamsToSubParams(const std::string name, const Teuchos::ParameterList& inparams)
    {
      (*params_).sublist(name) = TranslateSolverParameters(inparams);
    }
    //@}

    //! @name Query methods

    /*!
    \brief Get communicator

    */
    inline const Epetra_Comm& Comm() const { return comm_; }

    /*!
    \brief Get solver parameters

    */
    inline Teuchos::ParameterList& Params() const { return *params_; }

    //@}

    /*!
    \brief Return the solver name from the solver block in the input file

    \note This name is purely descriptive and does not affect any computations.
    */
    inline std::string Name() const { return params_->get<std::string>("name"); }

    /*!
    \brief Return number of iterations performed by solver
    */
    int getNumIters() const;

   private:
    /*!
   \brief Adapt tolerance of iterative solver

    This method allows to adapt the tolerance of the underlying iterative solver,
    if an iterative solver is used. It is meant to be used together with
    a relative convergence criteria AZ_r0 (decided from input file)
    and allows to adapt this relative convergence criteria depending on
    the current residual of the outer nonlinear solver

    It computes a new relative tolerance to be<br>
    <br>
    \code
    if (currentnlnres*tol < desirednlnres)
      tol = desirednlnres * better / currentnlnres
    \endcode

    \note This is a rule of thumb method - not a true adaptivity in the
          field of inexact Newton methods.

    \note This method works with iterative solvers only - it does nothing for all other
          solvers configured.

    \sa ResetTolerance

    \param desirednlnres (in): Desired residual in outer nonlinear solve
    \param currentnlnres (in): Current residual in outer nonlinear solve
    \param better        (in): The amount the linear solver shall be better than
                               currentnlnres

    */
    void AdaptTolerance(
        const double desirednlnres, const double currentnlnres, const double better);

    //! set tolerance to Belos solver
    void SetTolerance(double tolerance);

    //! a communicator
    const Epetra_Comm& comm_;
    //! (internal) parameter list
    Teuchos::RCP<Teuchos::ParameterList> params_;

    /// internal solver strategy
    Teuchos::RCP<CORE::LINEAR_SOLVER::SolverType<Epetra_Operator, Epetra_MultiVector>> solver_;

   private:
    //! don't want = operator
    Solver operator=(const Solver& old) = delete;

    //! don't want cctor
    Solver(const Solver& old) = delete;

  };  // class Solver

}  // namespace CORE::LINALG

BACI_NAMESPACE_CLOSE

#endif  // SOLVER_LINALG_SOLVER_H
