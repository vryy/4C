/*!----------------------------------------------------------------------
\file linalg_solver.H

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
#include "Amesos_Umfpack.h"
#include "Amesos_Lapack.h"

// Trilinos is configured with SuperLUDIST only in the parallel version
#ifdef PARALLEL
#include "Amesos_Superludist.h"
#endif

#include "drt_dserror.H"
#include "../drt_inpar/inpar_solver.H"

using namespace std;
using namespace Teuchos;


// forward declarations
// spooles stuff
#ifdef PARALLEL
#ifdef SPOOLES_PACKAGE
struct _FrontMtx;
struct _InpMtx;
struct _DenseMtx;
struct _Graph;
struct _IVL;
struct _ETree;
struct _SubMtxManager;
struct _ChvManager;
struct _Chv;
struct _IV;
struct _DV;
struct _DenseMtx;
struct _SolveMap;
#endif
#endif

// forward declarations
// this should not be here
namespace FLD
{
class VM3_Solver;
}

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

  \param matrix (in/out): system of equations
  \param x      (in/out): initial guess on input, solution on output
  \param b      (in)    : right hand side vector
  \param refactor (in)  : flag indicating whether system should be refactorized
  \param reset  (in)    : flag indicating whether all data from previous solves should
                           be recalculated including preconditioners
  */
  void Solve(RCP<Epetra_Operator>  matrix,
             RCP<Epetra_Vector>    x,
             RCP<Epetra_Vector>    b,
             bool refactor,
             bool reset = false);

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
  inline int Ncall() const { return ncall_; }

  /*!
  \brief Get the factored_ flag indicating whether system has been factorized before

  */
  inline bool IsFactored() const { return factored_; }
  
  /*!
  \brief //! Returns the condition number estimate for the current preconditioner, if one exists, returns -1.0 if no estimate.

  */
  inline double Condest() const { return aztec_->Condest(); }

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
  \brief Solve using aztec

  */
  virtual void Solve_aztec(const bool reset);

  //! a communicator
  const Epetra_Comm&        comm_;
  //! (internal) parameter list
  RCP<ParameterList>        params_;
  //! file to write output to (proc 0 only, can be NULL on input)
  FILE*                     outfile_;
  //! system of equations
  RCP<Epetra_Operator>      A_;
  //! modified system of equations
  RCP<Epetra_CrsMatrix>     Aplus_;
  //! Preconditioner (currently an ML, Ifpack, MLAPI or SIMPLER operator)
  RCP<Epetra_Operator>      P_;
  //! system of equations used for preconditioning used by P_ only
  RCP<Epetra_RowMatrix>     Pmatrix_;
  //! initial guess and solution
  RCP<Epetra_Vector>        x_;
  //! right hand side vector
  RCP<Epetra_Vector>        b_;
  //int x_LDA_;
  //int b_LDA_;
  //! a linear problem wrapper class used by Trilinos and for scaling of the system
  RCP<Epetra_LinearProblem> lp_;
  //! Epetra multivector: initial guess and solution
  Epetra_MultiVector *              X_;
  //! Epetra multivector: right hand side
  Epetra_MultiVector *              B_;
  //! flag indicating whether matrix was factored before
  bool                              factored_;
  //! counting how many times matrix was solved between resets
  int                               ncall_;
  //! an abstract amesos solver that can be any of the amesos concrete immplementations
  RCP<Amesos_BaseSolver>         amesos_;
  //! an aztec solver
  RCP<AztecOO>                    aztec_;
  //! an aztex maxiter status test
  RCP<AztecOO_StatusTestMaxIters> aztest_maxiter_;
  //! an aztec 2-norm convergence test
  RCP<AztecOO_StatusTestResNorm>  aztest_norm2_;
  //! an aztec inf-norm convergence test
  RCP<AztecOO_StatusTestResNorm>  aztest_norminf_;
  // ! combination of the above tests
  RCP<AztecOO_StatusTestCombo>    aztest_combo1_;
  RCP<AztecOO_StatusTestCombo>    aztest_combo2_;
  //! a vm3 solver
  RCP<FLD::VM3_Solver>            vm3_solver_;

  //! reindex linear problem for amesos
  RCP<EpetraExt::LinearProblem_Reindex> reindexer_;

  //! downwinding reindexer for downwind smoothers
  RCP<LINALG::DownwindMatrix>           dwind_;

#ifdef PARALLEL
#ifdef SPOOLES_PACKAGE
  //! lots of stuff needed by spooles
  _FrontMtx               *frontmtx_;
  _InpMtx                 *newA_;
  _DenseMtx               *newY_;
  _ETree                  *frontETree_;
  _SubMtxManager          *mtxmanager_;
  _IV                     *newToOldIV_;
  _IV                     *oldToNewIV_;
  _IV                     *ownersIV_;
  _IV                     *vtxmapIV_;
  _IV                     *ownedColumnsIV_;
  _SolveMap               *solvemap_;
  _IVL                    *symbfacIVL_;
  _Graph                  *graph_;
  _DV                     *cumopsDV_;
  _DenseMtx               *mtxY_;
  _ChvManager             *chvmanager_;
  _DenseMtx               *mtxX_;
  _InpMtx                 *mtxA_;
#endif
#endif


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

\param A (in)       : Matrix A for C=A*B
\param B (in)       : Matrix B for C=A*B
\param complete (in): flag indicating whether C shall be called FillComplete.
                      There intentionally is NO default value for this here.
*/
RCP<SparseMatrix> MLMultiply(const Epetra_CrsMatrix& A,
                             const Epetra_CrsMatrix& B,
                             bool complete);

/*!
\brief Matrix-Matrix Multiply C = A*B using ML

\sa MLMultiply(const Epetra_CrsMatrix& A,const Epetra_CrsMatrix& B,bool complete);

*/
RCP<SparseMatrix> MLMultiply(const SparseMatrix& A,
                             const SparseMatrix& B,
                             bool complete);

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
