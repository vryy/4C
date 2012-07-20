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
#ifndef LINALG_SOLVER_H
#define LINALG_SOLVER_H

#include <string>  // this can probably removed

#include <Ifpack.h> // for condest

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_solver.H"



#include "solver_solvertype.H"
#include "solver_preconditionertype.H"

using namespace std;

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
  Solver(Teuchos::RCP<Teuchos::ParameterList> params, const Epetra_Comm& comm, FILE* outfile);


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
  \brief Setup system of equations

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
  void Setup(
    Teuchos::RCP<Epetra_Operator>     matrix                     ,
    Teuchos::RCP<Epetra_MultiVector>       x                     ,
    Teuchos::RCP<Epetra_MultiVector>       b                     ,
    bool                     refactor                            ,
    bool                     reset   = false                     ,
    Teuchos::RCP<Epetra_MultiVector>  weighted_basis_mean_w= Teuchos::null,
    Teuchos::RCP<Epetra_MultiVector>  kernel_c= Teuchos::null             ,
    bool                     project = false                          );

  /// solve linear system after setup has been done
  void Solve();

  /// Apply inverse of the preconditioner
  /*!
    For krylov solver, only the preconditioner is called here. Direct solvers
    do a direct solve. Setup needs to be called before.
   */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  /*!
  \brief Solve system of equations in one go

  This includes setup. Reuse of preconditioners and factorized systems is
  provided.

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
      Teuchos::RCP<Epetra_Operator>     matrix                     ,
      Teuchos::RCP<Epetra_MultiVector>       x                     ,
      Teuchos::RCP<Epetra_MultiVector>       b                     ,
      bool                     refactor                            ,
      bool                     reset   = false                     ,
      Teuchos::RCP<Epetra_MultiVector>  weighted_basis_mean_w= Teuchos::null,
      Teuchos::RCP<Epetra_MultiVector>  kernel_c= Teuchos::null             ,
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
                        e.g. DRT::Problem::SolverParams(num) in case of solver for
                        structures and num according to STRUCTURAL DYNAMIC
  \return             : internal parameter list ready to be associated
                        with #params_
  \date 02/07,11/08
  */
  static const Teuchos::ParameterList TranslateSolverParameters(const Teuchos::ParameterList& inparams);

  /*!
   \brief Translate BACI dat file parameters to Stratimikos
   */
  static const Teuchos::ParameterList TranslateToStratimikos(const Teuchos::ParameterList& inparams);

  /*!
   \brief Translate BACI dat file parameters to Stratimikos (Amesos)
   */
  static const Teuchos::ParameterList TranslateAmesosToStratimikos(const Teuchos::ParameterList& inparams);

  static const Teuchos::ParameterList TranslateAztecToStratimikos(const Teuchos::ParameterList& inparams);

  static const Teuchos::ParameterList TranslateBelosToStratimikos(const Teuchos::ParameterList& inparams);

  static const Teuchos::ParameterList TranslateIfpackToStratimikos(const Teuchos::ParameterList& inparams);

  static const Teuchos::ParameterList TranslateMLToStratimikos(const Teuchos::ParameterList& inparams);

  static const Teuchos::ParameterList TranslateBACIToML(const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist = NULL);

  static const Teuchos::ParameterList TranslateBACIToIfpack(const Teuchos::ParameterList& inparams);

  static const Teuchos::ParameterList TranslateBACIToTeko(const Teuchos::ParameterList& inparams);

  static const Teuchos::ParameterList TranslateBACIToBelos(const Teuchos::ParameterList& inparams);

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

  /*!
  \brief Fix the nullspace of ML to match a new given map

  The nullspace is looked for in the parameter list. If found, it is assumed that
  it matches the oldmap. Then it is fixed to match the new map

  \param field (in): field name (just used for output)
  \param oldmap (in): row map of nullspace
  \param newmap (in): row map of nullspace upon exit
  \param solveparams (in): parameterlist including nullspace vector
  */
  static void FixMLNullspace(std::string field,
                      const Epetra_Map& oldmap,
                      const Epetra_Map& newmap,
                      Teuchos::ParameterList& solveparams);

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
  inline Teuchos::ParameterList& Params() const { return *params_; }

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

  /*!
  \brief //! Returns the solver name that has been set in the corresponding block in the BACI dat file
             used for debugging purposes only. has no meaning
  */
   inline std::string Name() const { return params_->get<std::string>("name"); }

private:
   //! @name BuildSolver methods

    /*!
    \brief helper function for generating BELOS solver

    creates an internal BelosSolver object and stores it in solver_
    */
    void BuildBelosSolver(const Epetra_Comm & comm, Teuchos::ParameterList & params, FILE * outfile);

  //@}

  //! a communicator
  const Epetra_Comm&                    comm_;
  //! (internal) parameter list
  Teuchos::RCP<Teuchos::ParameterList>                    params_;
  //! file to write output to (proc 0 only, can be NULL on input)
  FILE*                                 outfile_;

  /// internal solver strategy
  Teuchos::RCP<SOLVER::SolverType> solver_;

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
Teuchos::RCP<SparseMatrix> MLMultiply(const Epetra_CrsMatrix& Aorig,
                             const Epetra_CrsMatrix& Borig,
                             bool explicitdirichlet,
                             bool savegraph,
                             bool complete);

/*!
\brief Matrix-Matrix Multiply C = A*B using ML

\sa MLMultiply(const Epetra_CrsMatrix& A,const Epetra_CrsMatrix& B,bool complete);

*/
Teuchos::RCP<SparseMatrix> MLMultiply(const SparseMatrix& A,
                             const SparseMatrix& B,
                             bool explicitdirichlet,
                             bool savegraph,
                             bool complete);

/*!
\brief Matrix-Matrix Multiply C = A*B using ML

\sa MLMultiply(const Epetra_CrsMatrix& A,const Epetra_CrsMatrix& B,bool complete);

*/
Teuchos::RCP<SparseMatrix> MLMultiply(const SparseMatrix& A,
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
Teuchos::RCP<SparseMatrix> MLMultiply(const SparseMatrix& A, bool transA,
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
