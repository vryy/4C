/*!----------------------------------------------------------------------
\file linalg_solver.H

\class LINALG::Solver

\brief A general solver interface to Trilinos solvers and spooles

- The input of parameters needs to be completely reworked (hiwi job)

- This class should implement Epetra_Operator in the future

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef LINALG_SOLVER_H
#define LINALG_SOLVER_H

#include "vector"
#include "string"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"

#include "AztecOO.h"

#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"

#include "Amesos_Klu.h"
#include "Amesos_Umfpack.h"
#include "Amesos_Lapack.h"

#include "../drt_fluid/vm3_solver.H"

// Trilinos is configured with SuperLUDIST only in the parallel version
#ifdef PARALLEL
#include "Amesos_Superludist.h"
#endif

#include "drt_dserror.H"

using namespace std;
using namespace Teuchos;


// forward declarations
struct _SOLVAR;
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


/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{

/*!
\class Solver

\brief A general solver interface to Trilinos solvers and spooles

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
  Solver(RefCountPtr<ParameterList> params, const Epetra_Comm& comm, FILE* outfile);


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
  void Solve(RefCountPtr<Epetra_CrsMatrix> matrix,
             RefCountPtr<Epetra_Vector>    x,
             RefCountPtr<Epetra_Vector>    b,
             bool refactor,
             bool reset = false);

  /*!
  \brief Solve system of equations with an Operator
  
  Solve system of equations, where the matrix is not necessarily a true matrix
  but some kind of directional derivative operator.
  Note that this can only be used in combination with Aztec and without scaling.
  All other solvers will throw an error.

  \param matrix (in/out): system of equations
  \param x      (in/out): initial guess on input, solution on output
  \param b      (in)    : right hand side vector
  \param refactor (in)  : flag indicating whether system should be refactorized
  \param reset  (in)    : flag indicating whether all data from previous solves should
                           be recalculated including preconditioners
  */
  void Solve(RefCountPtr<Epetra_Operator> Operator,
             RefCountPtr<Epetra_Vector>    x,
             RefCountPtr<Epetra_Vector>    b,
             bool refactor,
             bool reset = false);

  /*!
  \brief Solve system of equations using VM3 solver (overloaded)

  \param matrix (in/out): system of equations
  \param x      (in/out): initial guess on input, solution on output
  \param b      (in)    : right hand side vector
  \param refactor (in)  : flag indicating whether system should be refactorized
  \param reset  (in)    : flag indicating whether all data from previous solves should
                           be recalculated including preconditioners
  \param matrix2 (in/out): mod. system of equations on finest grid level
  */
  void Solve(RefCountPtr<Epetra_CrsMatrix> matrix,
             RefCountPtr<Epetra_CrsMatrix> matrix2,
             RefCountPtr<Epetra_Vector>    x,
             RefCountPtr<Epetra_Vector>    b,
             bool refactor,
             bool reset = false);

  /*!
  \brief Reset the solver and clear data

  All data is destroyed except the parameter list
  */
  void Reset();

  //@}
  //! @name Input of parameters

  /*!
  \brief Translate solver input parameters from old ccarat style to
         new solver parameters list style
  
  This method should go away once the input of all solver parameters is reworked
         
  \param params (out): parameter list holding new style solver parameters
                       on output
  \param actsolv (in): ccarat style struct holding solver parameters

  */
  void TranslateSolverParameters(ParameterList& params, struct _SOLVAR* actsolv) const;

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
  inline const int Ncall() const { return ncall_; }

  /*!
  \brief Get the factored_ flag indicating whether system has been factorized before

  */
  inline const bool IsFactored() const { return factored_; }

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
  const Epetra_Comm&                comm_;    
  //! parameter list
  RefCountPtr<ParameterList>        params_;  
  //! file to write output to (proc 0 only, can be NULL on input)
  FILE*                             outfile_; 
  //! system of equations
  RefCountPtr<Epetra_Operator>      A_;       
  //! modified system of equations
  RefCountPtr<Epetra_CrsMatrix>     Aplus_;       
  //! Preconditioner (currently an ML or Ifpack operator)
  RefCountPtr<Epetra_Operator>      P_;       
  //! system of equations used for preconditioning used by P_ only
  RefCountPtr<Epetra_CrsMatrix>     Pmatrix_; 
  //! initial guess and solution
  RefCountPtr<Epetra_Vector>        x_;       
  //! right hand side vector
  RefCountPtr<Epetra_Vector>        b_;       
  //int x_LDA_;
  //int b_LDA_;
  //! a linear problem wrapper class used by Trilinos and for scaling of the system
  RefCountPtr<Epetra_LinearProblem> lp_;      
  //! Epetra multivector: initial guess and solution
  Epetra_MultiVector *              X_;
  //! Epetra multivector: right hand side
  Epetra_MultiVector *              B_;
  //! flag indicating whether matrix was factored before
  bool                              factored_;
  //! counting how many times matrix was solved between resets
  int                               ncall_;   
  //! an abstract amesos solver that can be any of the amesos concrete immplementations
  RefCountPtr<Amesos_BaseSolver>    amesos_;  
  //! an aztec solver
  RefCountPtr<AztecOO>              aztec_;  
  //! a vm3 solver
  RefCountPtr<VM3_Solver>           vm3_solver_; 

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
} // namespace LINALG


// << operator
ostream& operator << (ostream& os, const LINALG::Solver& solver);









#endif  // #ifndef LINALG_SOLVER_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
