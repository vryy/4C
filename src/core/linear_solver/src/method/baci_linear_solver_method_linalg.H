/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of basic solver interface

\level 0

*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_METHOD_LINALG_H
#define BACI_LINEAR_SOLVER_METHOD_LINALG_H

#include "baci_inpar_solver.H"
#include "baci_utils_exceptions.H"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_LinearProblem;

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
  /*!
  \brief A general solver interface to Trilinos solvers and spooles

  - The input of parameters needs to be completely reworked (hiwi job)

  - This class should implement Epetra_Operator in the future

  */
  class Solver
  {
   public:
    /*!
    \brief Standard Constructor

    \param params (in) : parameter list holding solver configuration
    \param comm (in)   : a reference to a epetra communicator object
    \param outfile (in): ptr to an error file on proc 0 only. If NULL on input,
                         no output will be generated.
    */
    Solver(Teuchos::RCP<Teuchos::ParameterList>& params, const Epetra_Comm& comm, FILE* outfile);

    /*!
    \brief Destructor

    */
    virtual ~Solver();

    /*!
    \brief Constructor with default solver Amesos_KLU

    Creates a default solver using nonsymmetric KLU without passing a
    parameter list. This is useful for quick development.
    default solver intentionally is not Umfpack as Umfpack is not a
     Trilinos build-in library.

    \param comm (in)   : a reference to a epetra communicator object
    \param outfile (in): ptr to an error file on proc 0 only. If NULL on input,
                         no output will be generated.
    */
    Solver(const Epetra_Comm& comm, FILE* outfile = NULL);

    /*!
    \brief Constructor taking a validated input parameter list for Solver

    Creates a solver using the parameters provided by inparams. This parameter
    list differs from the internally stored #params_. They are translated
    by #TranslateSolverParameters.

    \param inparams (in): input parameter list as provided by DRT::Problem,
                          e.g. DRT::Problem::SolverParams(num) in case of solver for
                          structures and num according to STRUCTURAL DYNAMIC
    \param comm     (in): a reference to a Epetra communicator object
    \param outfile  (in): pointer to an error file

    */
    Solver(const Teuchos::ParameterList& inparams, const Epetra_Comm& comm, FILE* outfile);

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
    \param refactor  (in)    : flag indicating whether system should be refactorized
    \param reset     (in)    : flag indicating whether all data from previous solves
                               should be recalculated including preconditioners
    \param projector (in)    : Krylov space projector - Teuchos::null is no
                               projection desired
    */
    void Setup(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
        Teuchos::RCP<Epetra_MultiVector> b, bool refactor, bool reset = false,
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector = Teuchos::null);

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
    \param refactor  (in)    : flag indicating whether system should be refactorized
    \param reset     (in)    : flag indicating whether all data from previous solves
                               should be recalculated including preconditioners
    \param projector (in)    : Krylov space projector - Teuchos::null is no
                               projection desired
    */
    int Solve(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
        Teuchos::RCP<Epetra_MultiVector> b, bool refactor, bool reset = false,
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector = Teuchos::null);

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
    \param refactor   (in)    : flag indicating whether system should be refactorized
    \param reset      (in)    : flag indicating whether all data from previous solves
                                should be recalculated including preconditioners
    \param projector  (in)    : Krylov space projector - Teuchos::null is no
                                projection desired
    */
    int NoxSolve(Epetra_LinearProblem& linProblem, bool refactor, bool reset = false,
        Teuchos::RCP<CORE::LINALG::KrylovProjector> projector = Teuchos::null);

    /*!
    \brief Reset the solver and clear data

    All data is destroyed except the parameter list
    */
    void Reset();

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

    \param inparams (in): input parameter list as provided by DRT::Problem,
                          e.g. DRT::Problem::SolverParams(num) in case of solver for
                          structures and num according to STRUCTURAL DYNAMIC
    \return             : internal parameter list ready to be associated
                          with #params_
    \date 02/07,11/08
    */
    static const Teuchos::ParameterList TranslateSolverParameters(
        const Teuchos::ParameterList& inparams);

    /*!
    \brief Translate BACI dat file parameters

    */
    static const Teuchos::ParameterList TranslateBACIToML(
        const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist = nullptr);

    static const Teuchos::ParameterList TranslateBACIToMuelu(
        const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist = nullptr);

    static const Teuchos::ParameterList TranslateBACIToIfpack(
        const Teuchos::ParameterList& inparams);

    static const Teuchos::ParameterList TranslateBACIToBelos(
        const Teuchos::ParameterList& inparams);

    /*!
    \brief Add a validated input parameter list as sublist to internal
           parameters

    \param name     (in): name for sublist for #params_ to attach data to
    \param inparams (in): input parameter list as provided by DRT::Problem,
                          e.g. DRT::Problem::FluidPressureSolverParams in case
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
    //! a communicator
    const Epetra_Comm& comm_;
    //! (internal) parameter list
    Teuchos::RCP<Teuchos::ParameterList> params_;
    //! file to write output to (proc 0 only, can be NULL on input)
    FILE* outfile_;

    /// internal solver strategy
    Teuchos::RCP<CORE::LINEAR_SOLVER::SolverType<Epetra_Operator, Epetra_MultiVector>> solver_;

   private:
    //! don't want = operator
    Solver operator=(const Solver& old) = delete;

    //! don't want cctor
    Solver(const Solver& old) = delete;

  };  // class Solver

}  // namespace CORE::LINALG

#endif  // SOLVER_LINALG_SOLVER_H
