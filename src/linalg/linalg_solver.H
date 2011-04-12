/*!----------------------------------------------------------------------
\file  linalg_solver.H

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef LINALG_SOLVER_H
#define LINALG_SOLVER_H

#include <vector>
#include <string>

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include <EpetraExt_Reindex_LinearProblem.h>

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"

#include "AztecOO.h"
#include "AztecOO_StatusTestResNorm.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestMaxIters.h"

#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"

#include "Amesos_Klu.h"
#ifndef HAVENOT_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#include "Amesos_Lapack.h"

// Trilinos is configured with SuperLUDIST only in the parallel version (not on altix)
#ifdef PARALLEL
#ifndef HAVENOT_SUPERLU
#include "Amesos_Superludist.h"
#endif
#endif

#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_solver.H"

#include "linalg_projected_precond.H"
#include "linalg_projected_operator.H"



using namespace std;
using namespace Teuchos;


namespace LINALG
{

// forward declarations
class SparseMatrix;
class SparseOperator;
class DownwindMatrix;

/*!
\brief A general solver interface to Trilinos solvers and spooles

- The input of parameters needs to be completely reworked (hiwi job)

- This class should implement Epetra_Operator in the future

\author gee (gee@lnm.mw.tum.de)
*/
class Solver
{

  /// preconditioner base class
  /*!

    A KrylovSolver object needs one (or more) preconditioner objects. There
    are many possible preconditioners. A unified framework simplifies the
    solution process.

    \author u.kue
    \date 04/11
   */
  class PreconditionerType
  {
  public:

    /// construction with output on proc 0 (might be NULL)
    /*!
      No setup is done upon construction, only the preconditioner object is
      created.
     */
    explicit PreconditionerType( FILE * outfile );

    /// virtual destruction
    virtual ~PreconditionerType() {}

    /// linear problem created (managed) by this preconditioner
    /*!
      This is how aztec sees the linear problem that needs to be solved.
     */
    virtual Epetra_LinearProblem & LinearProblem() { return *lp_; }

    /// support routine for setup
    void SetupLinearProblem( Epetra_Operator * matrix,
                             Epetra_MultiVector * x,
                             Epetra_MultiVector * b );

    /// Setup preconditioner with a given linear system.
    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b ) = 0;

    /// Finish calculation after linear solve.
    /*!
      This is empty in most cases, however some preconditioners might want to
      scale the solution.
     */
    virtual void Finish( Epetra_Operator * matrix,
                         Epetra_MultiVector * x,
                         Epetra_MultiVector * b ) {}

    /// linear operator used for preconditioning
    /*!
      This is how aztec sees the linear operator to be used for preconditioning.
     */
    virtual Epetra_Operator * PrecOperator() = 0;

  protected:

    //! a linear problem wrapper class used by Trilinos and for scaling of the system
    Teuchos::RCP<Epetra_LinearProblem> lp_;

    //! file to write output to (proc 0 only, can be NULL on input)
    FILE * outfile_;
  };

  /// ifpack preconditioners
  /*!
    Set of standard single-matrix preconditioners.
   */
  class IFPACKPreconditioner : public PreconditionerType
  {
  public:

    IFPACKPreconditioner( FILE * outfile, ParameterList & ifpacklist, ParameterList & azlist );

    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b );

    virtual Epetra_Operator * PrecOperator() { return &*prec_; }

  private:

    ParameterList & ifpacklist_;
    ParameterList & azlist_;

    //! system of equations used for preconditioning used by P_ only
    RCP<Epetra_RowMatrix>                 Pmatrix_;

    /// preconditioner
    RCP<Ifpack_Preconditioner> prec_;
  };

  /// ml preconditioners
  /*!
    Set of single-matrix algebraic multi-grid preconditioners.
   */
  class MLPreconditioner : public PreconditionerType
  {
  public:

    MLPreconditioner( FILE * outfile, ParameterList & mllist );

    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b );

    virtual Epetra_Operator * PrecOperator() { return &*P_; }

  private:

    ParameterList & mllist_;

    //! system of equations used for preconditioning used by P_ only
    RCP<Epetra_RowMatrix>                 Pmatrix_;

    /// preconditioner
    RCP<Epetra_Operator>                  P_;

  };

  /// downwinding, a preconditioner-wrapper
  /*!
    Downwinding does a matrix-reordering that is supposed to speed-up the
    solution process. The actual solution is performed using another
    preconditioner. DWindPreconditioner is just a wrapper that presents a
    permuted matrix to the actual preconditioner (and solver).

    \note DWindPreconditioner is a preconditioner in the solver framework, not
    in a mathematical sense.
   */
  class DWindPreconditioner : public PreconditionerType
  {
  public:

    DWindPreconditioner( FILE * outfile,
                         Teuchos::RCP<PreconditionerType> preconditioner,
                         ParameterList & azlist );

    /// return the linear problem of the wrapped preconditioner, thus showing
    /// the permuted system to aztec.
    virtual Epetra_LinearProblem & LinearProblem()
    {
      return preconditioner_->LinearProblem();
    }

    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b );

    /// undo reordering of lhs, don't care for rhs or matrix
    virtual void Finish( Epetra_Operator * matrix,
                         Epetra_MultiVector * x,
                         Epetra_MultiVector * b );

    /// return the wrapped preconditioner, no work needs to be done after the
    /// permutation is performed
    virtual Epetra_Operator * PrecOperator() { return preconditioner_->PrecOperator(); }

  private:

    Teuchos::RCP<PreconditionerType> preconditioner_;

    ParameterList & azlist_;

    //! downwinding reindexer for downwind smoothers
    RCP<LINALG::DownwindMatrix>            dwind_;

    RCP<Epetra_CrsMatrix>     dwA_;
    RCP<Epetra_MultiVector>   dwx_;
    RCP<Epetra_MultiVector>   dwb_;
  };

  /// SIMPLE(R) block preconditioner
  /*!
    Block preconditioners assume the Epetra_Operator to be a
    LINALG::BlockSparseMatrix.
   */
  class SimplePreconditioner : public PreconditionerType
  {
  public:

    SimplePreconditioner( FILE * outfile,
                          ParameterList & params, ParameterList & simpleparams );

    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b );

    virtual Epetra_Operator * PrecOperator() { return &*P_; }

  private:
    ParameterList & params_;
    ParameterList & simpleparams_;

    RCP<Epetra_Operator>    P_;
  };

  /// AMG BS
  class AMGBSPreconditioner : public PreconditionerType
  {
  public:

    AMGBSPreconditioner( FILE * outfile,
                         ParameterList & params );

    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b );

    virtual Epetra_Operator * PrecOperator() { return &*P_; }

  private:
    ParameterList & params_;
    RCP<Epetra_Operator>    P_;
  };

  /// General purpose block gauss-seidel preconditioner
  /*!
    2x2 block preconditioner
   */
  class BGSPreconditioner : public PreconditionerType
  {
  public:

    BGSPreconditioner( FILE * outfile,
                       ParameterList & params, ParameterList & bgslist );

    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b );

    virtual Epetra_Operator * PrecOperator() { return &*P_; }

  private:

    ParameterList & params_;
    ParameterList & bgslist_;
    RCP<Epetra_Operator>    P_;
  };

  /// krylov projection for undefined pressure value in incompressible fluids
  /*!
    This is not a preconditioner in a mathematical sense, but it fits the
    software framework nicely.

    A "real" preconditioner is wrapped.
   */
  class KrylovProjectionPreconditioner : public PreconditionerType
  {
  public:

    KrylovProjectionPreconditioner( FILE * outfile,
                                    Teuchos::RCP<PreconditionerType> preconditioner,
                                    RCP<Epetra_MultiVector>  weighted_basis_mean,
                                    RCP<Epetra_MultiVector>  kernel_c );

    virtual Epetra_LinearProblem & LinearProblem() { return preconditioner_->LinearProblem(); }

    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b );

    virtual void Finish( Epetra_Operator * matrix,
                         Epetra_MultiVector * x,
                         Epetra_MultiVector * b );

    /// return the projecting operator
    virtual Epetra_Operator * PrecOperator() { return &*A_; }

  private:
    Teuchos::RCP<PreconditionerType> preconditioner_;

    /// Peter's projector object that does the actual work
    Teuchos::RCP<LINALG::KrylovProjector> projector_;
    Teuchos::RCP<Epetra_MultiVector> weighted_basis_mean_;
    Teuchos::RCP<Epetra_MultiVector> kernel_c_;

    /// linear operator that calls a "real" preconditioning operator and does
    /// a projection afterwards.
    Teuchos::RCP<LINALG::LinalgProjectedOperator> A_;
  };

  /// inf-norm scaling fits the preconditioner framework perfecty
  /*!
    Modifies the underlying matrix and needs to unscale the result
    afterwards. Can be combined with any single-matrix preconditioner.
   */
  class InfNormPreconditioner : public PreconditionerType
  {
  public:

    InfNormPreconditioner( Teuchos::RCP<PreconditionerType> preconditioner );

    virtual Epetra_LinearProblem & LinearProblem() { return preconditioner_->LinearProblem(); }

    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b );

    virtual void Finish( Epetra_Operator * matrix,
                         Epetra_MultiVector * x,
                         Epetra_MultiVector * b );

    virtual Epetra_Operator * PrecOperator() { return preconditioner_->PrecOperator(); }

  private:
    Teuchos::RCP<PreconditionerType> preconditioner_;
    RCP<Epetra_Vector> rowsum_;
    RCP<Epetra_Vector> colsum_;
  };

  /// diagonal scaling fits the preconditioner framework perfecty
  /*!
    Modifies the underlying matrix and needs to unscale the result
    afterwards. Can be combined with any single-matrix preconditioner.
   */
  class SymDiagPreconditioner : public PreconditionerType
  {
  public:

    SymDiagPreconditioner( Teuchos::RCP<PreconditionerType> preconditioner );

    virtual Epetra_LinearProblem & LinearProblem() { return preconditioner_->LinearProblem(); }

    virtual void Setup( bool create,
                        Epetra_Operator * matrix,
                        Epetra_MultiVector * x,
                        Epetra_MultiVector * b );

    virtual void Finish( Epetra_Operator * matrix,
                         Epetra_MultiVector * x,
                         Epetra_MultiVector * b );

    virtual Epetra_Operator * PrecOperator() { return preconditioner_->PrecOperator(); }

  private:

    /// embedded preconditioner
    Teuchos::RCP<PreconditionerType> preconditioner_;

    /// matrix diagonal
    RCP<Epetra_Vector> diag_;
  };

  /// linear solver type base class
  class SolverType
  {
  public:

    virtual ~SolverType() {}

    /// solve the linear system
    virtual void Solve( RCP<Epetra_Operator>     matrix             ,
                        RCP<Epetra_Vector>       x                  ,
                        RCP<Epetra_Vector>       b                  ,
                        bool                             refactor           ,
                        bool                             reset              ,
                        RCP<Epetra_MultiVector>  weighted_basis_mean,
                        RCP<Epetra_MultiVector>  kernel_c           ,
                        bool                             project) = 0;

    /// preconditioner used within this solver
    virtual PreconditionerType * Preconditioner() = 0;
  };

  /// direct linear solver (using amesos)
  class DirectSolver : public SolverType
  {
  public:

    explicit DirectSolver( std::string solvertype );

    ~DirectSolver();

    virtual void Solve( RCP<Epetra_Operator>     matrix             ,
                        RCP<Epetra_Vector>       x                  ,
                        RCP<Epetra_Vector>       b                  ,
                        bool                             refactor           ,
                        bool                             reset              ,
                        RCP<Epetra_MultiVector>  weighted_basis_mean,
                        RCP<Epetra_MultiVector>  kernel_c           ,
                        bool                             project);

    virtual PreconditionerType * Preconditioner() { return NULL; }

    bool IsFactored() { return factored_; }

  private:

    std::string solvertype_;

    //! flag indicating whether matrix was factored before
    bool                                  factored_;

    //! a linear problem wrapper class used by Trilinos and for scaling of the system
    RCP<Epetra_LinearProblem>             lp_;

    //! initial guess and solution
    RCP<Epetra_Vector>                    x_;

    //! right hand side vector
    RCP<Epetra_Vector>                    b_;

    //! system of equations
    RCP<Epetra_Operator>  A_;

    //! an abstract amesos solver that can be any of the amesos concrete implementations
    RCP<Amesos_BaseSolver>                amesos_;

    //! reindex linear problem for amesos
    RCP<EpetraExt::LinearProblem_Reindex> reindexer_;
  };

  /// kryolv subspace linear solvers (aztec) with right-side preconditioning
  class KrylovSolver : public SolverType
  {
  public:

    KrylovSolver( const Epetra_Comm & comm, Teuchos::ParameterList & params, FILE * outfile );

    ~KrylovSolver();

    virtual void Solve( RCP<Epetra_Operator>     matrix             ,
                        RCP<Epetra_Vector>       x                  ,
                        RCP<Epetra_Vector>       b                  ,
                        bool                             refactor           ,
                        bool                             reset              ,
                        RCP<Epetra_MultiVector>  weighted_basis_mean,
                        RCP<Epetra_MultiVector>  kernel_c           ,
                        bool                             project);

    virtual PreconditionerType * Preconditioner() { return &*preconditioner_; }

    int Ncall() { return ncall_; }

    ParameterList & Params() const { return params_; }

  private:

    void CreatePreconditioner( ParameterList & azlist,
                               bool isCrsMatrix,
                               RCP<Epetra_MultiVector> weighted_basis_mean,
                               RCP<Epetra_MultiVector> kernel_c,
                               bool project );

    //! a communicator
    const Epetra_Comm&                    comm_;

    //! (internal) parameter list
    ParameterList&                  params_;

    //! file to write output to (proc 0 only, can be NULL on input)
    FILE*                                 outfile_;

    //! initial guess and solution
    RCP<Epetra_Vector>                    x_;

    //! right hand side vector
    RCP<Epetra_Vector>                    b_;

    //! system of equations
    RCP<Epetra_Operator>  A_;

    //! counting how many times matrix was solved between resets
    int                                   ncall_;

    Teuchos::RCP<PreconditionerType> preconditioner_;

    //! an aztex maxiter status test
    RCP<AztecOO_StatusTestMaxIters>       aztest_maxiter_;
    //! an aztec 2-norm convergence test
    RCP<AztecOO_StatusTestResNorm>        aztest_norm2_;
    //! an aztec inf-norm convergence test
    RCP<AztecOO_StatusTestResNorm>        aztest_norminf_;
    // ! combination of the above tests
    RCP<AztecOO_StatusTestCombo>          aztest_combo1_;
    RCP<AztecOO_StatusTestCombo>          aztest_combo2_;
  };

public:

  /*!
  \brief Standard Constructor

  \param params (in) : parameter list holding solver configuration
  \param comm (in)   : a reference to a epetra communicator object
  \param outfile (in): ptr to an error file on proc 0 only. If NULL on input,
                       no output will be generated.
  */
  Solver(RCP<ParameterList> params, const Epetra_Comm& comm, FILE* outfile);


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
                        e.g. DRT::Problem::StructSolverParams in case of solver for structures
  \param comm     (in): a reference to a Epetra communicator object
  \param outfile  (in): pointer to an error file

  \date 11/08
  */
  Solver(const Teuchos::ParameterList& inparams, const Epetra_Comm& comm, FILE* outfile);

  /*!
  \brief Set-up of stuff common to all constructors

  \date 11/08
  */
  void Setup();

  /*!
  \brief Destructor

  */
  virtual ~Solver();


  //! @name Solve and ApplyInverse methods

  /*!
  \brief Solve system of equations

  \param matrix                (in/out): system of equations
  \param x                     (in/out): initial guess on input, solution on output
  \param b                     (in/out): right hand side vector
                                         if project!=false it might be modified if not orthognal
                                         to matrix kernel.
  \param refactor              (in)    : flag indicating whether system should be refactorized
  \param reset                 (in)    : flag indicating whether all data from previous solves
                                         should be recalculated including preconditioners
  \param weighted_basis_mean_w (in)    : weight vector defining projector for iterative solution
                                         of singluar problems using iterative solvers. Required
                                         if project!=false. (Aztec only option)
                                         x will be orthogonalized to these weights on output.
  \param kernel_c              (in)    : Vector of ones in the respective components. Necessary
                                         if project!=false. (Aztec only option)
                                         Krylov space will be orthogonalized to this basis!
  \param project               (in)    : flag indicating whether several modes will be projected
                                         out of the Krylov space (Aztec only option)
  */
  void Solve(
    RCP<Epetra_Operator>     matrix                              ,
    RCP<Epetra_Vector>       x                                   ,
    RCP<Epetra_Vector>       b                                   ,
    bool                     refactor                            ,
    bool                     reset   = false                     ,
    RCP<Epetra_MultiVector>  weighted_basis_mean_w= Teuchos::null,
    RCP<Epetra_MultiVector>  kernel_c= Teuchos::null             ,
    bool                     project = false                     );

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

  \note This method works with Aztec only - it does nothing for all other
        solvers configured.

  \sa ResetTolerance

  \param desirednlnres (in): Desired residual in outer nonlinear solve
  \param currentnlnres (in): Current residual in outer nonlinear solve
  \param better        (in): The amount the linear solver shall be better than
                             currentnlnres

  */
  void AdaptTolerance(const double desirednlnres,
                      const double currentnlnres,
                      const double better);

  /*!
  \brief Adapt tolerance of iterative solver

  Reset the tolerance read from input file. Can only be used after a call to
  AdaptTolerance.

  \note This method works with Aztec only - it does nothing for all other
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
                        e.g. DRT::Problem::StructSolverParams in case of
                        solver for structures
  \return             : internal parameter list ready to be associated
                        with #params_
  \date 02/07,11/08
  */
  const Teuchos::ParameterList TranslateSolverParameters(const Teuchos::ParameterList& inparams);

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
  void PutSolverParamsToSubParams(const std::string name,
                                  const Teuchos::ParameterList& inparams)
  {
    (*params_).sublist(name) = TranslateSolverParameters(inparams);
    return;
  }
  //@}

  //! @name Query methods

  /*!
  \brief Print solver configuration (ostream << also supported)

  */
  virtual void Print(ostream& os) const;

  /*!
  \brief Get communicator

  */
  inline const Epetra_Comm& Comm() const { return comm_; }

  /*!
  \brief Get solver parameters

  */
  inline ParameterList& Params() const { return *params_; }

  /*!
  \brief Get number of solver calls done on this solver

  */
//   inline int Ncall() const { return ncall_; }

  /*!
  \brief Get the factored_ flag indicating whether system has been factorized before

  */
//   inline bool IsFactored() const { return factored_; }

  /*!
  \brief //! Returns the condition number estimate for the current preconditioner, if one exists, returns -1.0 if no estimate.

  */
//   inline double Condest() const { return aztec_->Condest(); }

  //@}

private:

  /*!
  \brief Solve using spooles

  */
  virtual void Solve_spooles(const bool reset);

  /*!
  \brief Solve using superlu

  */
  virtual void Solve_superlu(const bool reset);

  /*!
  \brief Solve using klu

  */
  virtual void Solve_klu(const bool reset);

  /*!
  \brief Solve using umfpack

  */
  virtual void Solve_umfpack(const bool reset);

  /*!
  \brief Solve using lapack

  */
  virtual void Solve_lapack(const bool reset);

  /*!
  \brief Solve using aztec. If project==false, system-matrix has to be
         non-singular. For project==true, the user can specify the
         matrix kernel and a weight vector to solve the problem
         iteratively on a residue class modulo the kernel.
         The representative returned in x_ will be orthogonal to
         the weighted_basis_mean_w vector.

  \param (i) reset
             flag indicating whether all data from previous solves
             should be recalculated including preconditioners
  \param (i) project
             flag indicating whether several modes will be projected
             out of the Krylov space (Aztec only option)
  \param (i) weighted_basis_mean_w
             weight vector defining projector for iterative solution
             of singluar problems using iterative solvers. Required
             if project!=false.
             x will be orthogonalized to these weights on output.
  \param (i) kernel_c
             Vector of ones in the respective components. Necess
             if project!=false.
             Krylov space will be orthogonalized to this basis!

  */
  virtual void Solve_aztec(
    const bool                         reset    ,
    const bool                         project  ,
    const RCP<LINALG::KrylovProjector> projector
    );

  //! a communicator
  const Epetra_Comm&                    comm_;
  //! (internal) parameter list
  RCP<ParameterList>                    params_;
  //! file to write output to (proc 0 only, can be NULL on input)
  FILE*                                 outfile_;

  RCP<SolverType> solver_;

private:

  // don't want = operator
  Solver operator = (const Solver& old);
  // don't want cctor
  Solver(const LINALG::Solver& old);

}; // class Solver


/*!
\brief Matrix-Matrix Multiply C = A*B using ML

Depending on the structure of your matrix, this method is potentially
(significantly) faster than the LINALG::Multiply implementations.
A factor of 5 - 10 in speed has been observed when comparing.
It uses ML instead of EpetraExt for the multiplication kernel.

<pre>
                       ***Warning***
         This method is not well tested (yet) and there is some risk of
         failure AND wrong results without notice. This depends
         on structural properties of the input matrices A and B.
         First:
         A,B should NOT contain columns in the
         column and/or domain map on a processor,
         where there is no nonzero entry in that column on a specific processor.
         This will definitely lead to failure but is a very rare case.
         Second:
         There should not be a row in row/rangemap with no
         nonzero entry in that row on any processor.
         Third:
         Very weirdo column orderings might lead to trouble, but
         scattered row/column ordering as usually the case in PDE
         problems should work fine.
         Fourth:
         If you have a choice how to formulate your problem:
</pre>

\note The Mat-Mat-Mult kernel used herein
      was optimized for algebraic multigrid. This means, the method
      will be significantly faster when B is sparser and/or smaller than A.

\note There is no natural transpose multiply in ML. So if you need a
      transpose-Multiply you have to transpose outside (or write a wrapper
      for this method)

\note ML wipes all exact zero entries in the product. So its perfectly
      ok if your product has less nonzeros than with the other
      LINALG::Multiply methods.

\param Aorig (in)       : Matrix A for C=A*B
\param Borig (in)       : Matrix B for C=A*B
\param complete (in): flag indicating whether C shall be called FillComplete.
                      There intentionally is NO default value for this here.
*/
RCP<SparseMatrix> MLMultiply(const Epetra_CrsMatrix& Aorig,
                             const Epetra_CrsMatrix& Borig,
                             bool explicitdirichlet,
                             bool savegraph,
                             bool complete);

/*!
\brief Matrix-Matrix Multiply C = A*B using ML

\sa MLMultiply(const Epetra_CrsMatrix& A,const Epetra_CrsMatrix& B,bool complete);

*/
RCP<SparseMatrix> MLMultiply(const SparseMatrix& A,
                             const SparseMatrix& B,
                             bool explicitdirichlet,
                             bool savegraph,
                             bool complete);

/*!
\brief Matrix-Matrix Multiply C = A*B using ML

\sa MLMultiply(const Epetra_CrsMatrix& A,const Epetra_CrsMatrix& B,bool complete);

*/
RCP<SparseMatrix> MLMultiply(const SparseMatrix& A,
                             const SparseMatrix& B,
                             bool complete);


/// Multiply a (transposed) matrix with another (transposed): C = A(^T)*B(^T)
/*!
  Multiply one matrix with another. Both matrices must be completed.
  Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

  \note This is a true parallel multiplication, even in the transposed case.

  \note Does call complete on C upon exit by default.

  \note Uses ML as multiplication kernel, not EpetraExt.

  \note In this version the flags explicitdirichlet and savegraph must be handed in.
        Thus, they can be defined explicitly, while in the standard version of Multipliy()
        above, result matrix C automatically inherits these flags from input matrix A

  \param A              (in)     : Matrix to multiply with B (must have Filled()==true)
  \param transA         (in)     : flag indicating whether transposed of A should be used
  \param B              (in)     : Matrix to multiply with A (must have Filled()==true)
  \param transB         (in)     : flag indicating whether transposed of B should be used
  \param explicitdirichlet (in)  : flag deciding on explicitdirichlet flag of C
  \param savegraph      (in)     : flag deciding on savegraph flag of C
  \param completeoutput (in)     : flag indicating whether Complete(...) shall be called on C upon output
  \return Matrix product A(^T)*B(^T)
*/
RCP<SparseMatrix> MLMultiply(const SparseMatrix& A, bool transA,
                             const SparseMatrix& B, bool transB,
                             bool explicitdirichlet,
                             bool savegraph,
                             bool completeoutput);

/*!
  \brief Compute a condition number estimate of given matrix A

  condition number kappa of a matrix A is defined as kappa(A) = ||A|| * ||A^-1||
  for A = singulaer, kappa = infinity

  \return kappa

  copied and adapted from
  Trilinos/packages/ifpack/src/Ifpack_Condest.h and Ifpack_Condest.cpp
*/
double Condest(
    SparseMatrix&                            A,  ///< sparse matrix
    const Ifpack_CondestType CT = Ifpack_Cheap,  ///< select approach for computing condition number estimate
    const int MaxIters = 1000,                   ///< maximum number of linear iterations
    const double Tol = 1e-9                      ///< tolerance
    );

} // namespace LINALG


// << operator
ostream& operator << (ostream& os, const LINALG::Solver& solver);









#endif  // #ifndef LINALG_SOLVER_H
#endif  // #ifdef CCADISCRET
