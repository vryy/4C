/*!----------------------------------------------------------------------
\file solver_cheapsimplepreconditioner.cpp

<pre>
Maintainer: Tobias Wiesner
            wiesner@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

// Trilinos headers
#include <EpetraExt_OperatorOut.h>
#include <Ifpack.h>
#include <ml_MultiLevelPreconditioner.h>

// BACI headers
#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"    // helper functions (linear Algebra related)
#include "../linalg/linalg_downwindmatrix.H"

#include "solver_cheapsimplepreconditioner.H"

#define SIMPLEC_DIAGONAL      1    // 1: row sums     0: just diagonal
#define CHEAPSIMPLE_ALGORITHM 1    // 1: AMG          0: true solve
#define SIMPLER_ALGORITHM     0    // 1: triple solve 0: double solve
#define SIMPLER_ALPHA         0.8  // simple pressure damping parameter
#define SIMPLER_TIMING        0    // printout timing of setup
/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 02/08|
 *----------------------------------------------------------------------*/
LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::CheapSIMPLE_BlockPreconditioner(RCP<Epetra_Operator> A,
                                           const ParameterList& predict_list,
                                           const ParameterList& correct_list,
                                           FILE* outfile)
  : outfile_(outfile),
    predictSolver_list_(predict_list),
    schurSolver_list_(correct_list),
    alpha_(SIMPLER_ALPHA),
    vdw_(false),
    pdw_(false),
    label_(SetupLabel())
{
  // remove the SIMPLER sublist from the predictSolver_list_,
  // otherwise it will try to recursively create a SIMPLE
  // preconditioner when we do the subblock solvers
  //if (predictSolver_list_.isSublist("SIMPLER")) predictSolver_list_.remove("SIMPLER");

  // check for contact, meshtying or constraints
  // (no special functionality yet, only checking)
//  const int myrank = A->Comm().MyPID();
//  bool mt = schurSolver_list_.get<bool>("MESHTYING",false);
//  if (!myrank && mt) cout << "\n**********\nMESHTYING SIMPLER\n**********\n\n";
//  bool co = schurSolver_list_.get<bool>("CONTACT",false);
//  if (!myrank && co) cout << "\n**********\nCONTACT SIMPLER\n**********\n\n";
//  bool cstr = schurSolver_list_.get<bool>("CONSTRAINT",false);
//  if (!myrank && cstr) cout << "\n**********\nCONSTRAINT SIMPLER\n**********\n\n";

  Setup(A,predict_list,correct_list);

  return;
}


/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::Setup(RCP<Epetra_Operator> A,
                                     const ParameterList& origvlist,
                                     const ParameterList& origplist)
{
  const int myrank = A->Comm().MyPID();
  Epetra_Time time(A->Comm());
  Epetra_Time totaltime(A->Comm());
  const bool visml = predictSolver_list_.isSublist("ML Parameters");
  const bool pisml = schurSolver_list_.isSublist("ML Parameters");
  const bool visifpack = predictSolver_list_.isSublist("IFPACK Parameters");
  const bool pisifpack = schurSolver_list_.isSublist("IFPACK Parameters");
  if (!visml && !visifpack) dserror("Have to use either ML or Ifpack for velocities");
  if (!pisml && !pisifpack) dserror("Have to use either ML or Ifpack for pressure");

  //-------------------------------------------------------------------------
  // either do manual split or use provided BlockSparseMatrixBase
  //-------------------------------------------------------------------------
  A_ = rcp_dynamic_cast<BlockSparseMatrixBase>(A);
  if (A_!=null)
  {
    // Make a shallow copy of the block matrix as the preconditioners on the
    // blocks will be reused and the next assembly will replace the block
    // matrices.
    A_ = A_->Clone(View);
    mmex_ = A_->RangeExtractor();
  }

  // if the MESHTYING, CONTACT or CONSTRAINT flag is set,
  // adapt multigrid nullspace for constraint block

  //-------------------------------------------------------------------------
  // Modify lists to reuse subblock preconditioner at least maxiter times
  //-------------------------------------------------------------------------
  {
    int maxiter = predictSolver_list_.sublist("Aztec Parameters").get("AZ_max_iter",1);
    predictSolver_list_.sublist("Aztec Parameters").set("reuse",maxiter+1);
  }

#if SIMPLEC_DIAGONAL
  //-------------------------------------------------------------------------
  // Allocate and compute abs(rowsum(A(0,0))^{-1}
  //-------------------------------------------------------------------------
  {
    Epetra_Vector diag(*mmex_.Map(0),false);
    RCP<Epetra_CrsMatrix> A00 = (*A_)(0,0).EpetraMatrix();
    A00->InvRowSums(diag);
    diagAinv_ = rcp(new SparseMatrix(diag));
    diagAinv_->Complete(*mmex_.Map(0),*mmex_.Map(0));
  }
#else
  //-------------------------------------------------------------------------
  // Allocate and compute diag(A(0,0)^{-1}
  //-------------------------------------------------------------------------
  {
    Epetra_Vector diag(*mmex_.Map(0),false);
    (*A_)(0,0).ExtractDiagonalCopy(diag);
    int err = diag.Reciprocal(diag);
    if (err) dserror("Epetra_MultiVector::Reciprocal returned %d",err);
    diagAinv_ = rcp(new SparseMatrix(diag));
    diagAinv_->Complete(*mmex_.Map(0),*mmex_.Map(0));
  }
#endif
  if (!myrank && SIMPLER_TIMING) printf("--- Time to do diagF^{-1}   %10.3E\n",time.ElapsedTime());
  time.ResetStartTime();

  //-------------------------------------------------------------------------
  // Allocate and compute approximate Schur complement operator S
  // S = A(1,1) - A(1,0) * diagAinv * A(0,1)
  //-------------------------------------------------------------------------
  {
    Epetra_Time ltime(A_->Comm());
    S_ = LINALG::Multiply(*diagAinv_,false,(*A_)(0,1),false,true);
    if (!myrank && SIMPLER_TIMING) printf("*** S = diagAinv * A(0,1) %10.3E\n",ltime.ElapsedTime());
    ltime.ResetStartTime();
    S_ = LINALG::MLMultiply((*A_)(1,0),*S_,false);
    // The LINALG::Multiply method would consume a HUGE amount of memory!!!
    // So always use LINALG::MLMultiply in the line above! Otherwise you won't be able
    // to solve any large linear problem since you'll definitely run out of memory.
    //S_ = LINALG::Multiply((*A_)(1,0),false,*S_,false,false);
    if (!myrank && SIMPLER_TIMING) printf("*** S = A(1,0) * S (ML)   %10.3E\n",ltime.ElapsedTime());
    ltime.ResetStartTime();
    S_->Add((*A_)(1,1),false,1.0,-1.0);
    if (!myrank && SIMPLER_TIMING) printf("*** S = A(1,1) - S        %10.3E\n",ltime.ElapsedTime());
    ltime.ResetStartTime();
    S_->Complete((*A_)(1,1).DomainMap(),(*A_)(1,1).RangeMap());
    if (!myrank && SIMPLER_TIMING) printf("*** S complete            %10.3E\n",ltime.ElapsedTime());
    ltime.ResetStartTime();


  }
  if (!myrank && SIMPLER_TIMING) printf("--- Time to do S            %10.3E\n",time.ElapsedTime());
  time.ResetStartTime();


#if CHEAPSIMPLE_ALGORITHM
  {
    Epetra_CrsMatrix* A00 = NULL;
    Epetra_CrsMatrix* A11 = NULL;
    A00 = (*A_)(0,0).EpetraMatrix().get();
    A11 = S_->EpetraMatrix().get();

  //-------------------------------------------------------------------------
  // Allocate preconditioner for pressure and velocity
  //-------------------------------------------------------------------------
    if (visml)
      Ppredict_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*A00,predictSolver_list_.sublist("ML Parameters"),true));
    else
    {
      string type = predictSolver_list_.sublist("Aztec Parameters").get("Preconditioner Type","ILU");
      Ifpack factory;
      Ifpack_Preconditioner* prec = factory.Create(type,A00,0);
      prec->SetParameters(predictSolver_list_.sublist("IFPACK Parameters"));
      prec->Initialize();
      prec->Compute();
      Ppredict_ = rcp(prec);
    }
    if (!myrank && SIMPLER_TIMING) printf("--- Time to do P(v)         %10.3E\n",time.ElapsedTime());
    time.ResetStartTime();

    if (pisml)
    {
      Pschur_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*A11,schurSolver_list_.sublist("ML Parameters"),true));
    }
    else
    {
      Ifpack factory;
      string type = schurSolver_list_.sublist("Aztec Parameters").get("Preconditioner Type","ILU");
      Ifpack_Preconditioner* prec = factory.Create(type,A11,0);
      prec->SetParameters(schurSolver_list_.sublist("IFPACK Parameters"));
      prec->Initialize();
      prec->Compute();
      Pschur_ = rcp(prec);
    }
    if (!myrank && SIMPLER_TIMING) printf("--- Time to do P(p)         %10.3E\n",time.ElapsedTime());
    time.ResetStartTime();
  }
#else
  //-------------------------------------------------------------------------
  // Allocate solver for pressure and velocity
  //-------------------------------------------------------------------------
  {
    RCP<ParameterList> vrcplist = rcp(&predictSolver_list_,false);
    vsolver_ = rcp(new LINALG::Solver(vrcplist,A_->Comm(),outfile_));
    RCP<ParameterList> prcplist = rcp(&schurSolver_list_,false);
    psolver_ = rcp(new LINALG::Solver(prcplist,A_->Comm(),outfile_));
  }
#endif

  //-------------------------------------------------------------------------
  // Allocate velocity and pressure solution and rhs vectors
  //-------------------------------------------------------------------------
  vx_ = rcp(new LINALG::ANA::Vector(*mmex_.Map(0),false));
  vb_ = rcp(new LINALG::ANA::Vector(*mmex_.Map(0),false));
  px_ = rcp(new LINALG::ANA::Vector(*mmex_.Map(1),false));
  pb_ = rcp(new LINALG::ANA::Vector(*mmex_.Map(1),false));

  //-------------------------------------------------------------------------
  // Allocate working vectors for velocity and pressure
  //-------------------------------------------------------------------------
  vwork1_  = rcp(new LINALG::ANA::Vector(*mmex_.Map(0),false));
  vwork2_  = rcp(new LINALG::ANA::Vector(*mmex_.Map(0),false));
  pwork1_ = rcp(new LINALG::ANA::Vector(*mmex_.Map(1),false));
  pwork2_ = rcp(new LINALG::ANA::Vector(*mmex_.Map(1),false));

  if (!myrank && SIMPLER_TIMING) printf("--- Time to do allocate mem %10.3E\n",time.ElapsedTime());
  if (!myrank && SIMPLER_TIMING) printf("=== Total simpler setup === %10.3E\n",totaltime.ElapsedTime());

  return;
}


/*----------------------------------------------------------------------*
 |  apply const operator (public)                            mwgee 02/08|
 *----------------------------------------------------------------------*/
int LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::ApplyInverse(const Epetra_MultiVector& X,
                                           Epetra_MultiVector& Y) const
{
  // note: Aztec might pass X and Y as physically identical objects,
  // so we better deep copy here

  // extract initial guess and rhs for velocity and pressure
  mmex_.ExtractVector(X,0,*vb_);
  mmex_.ExtractVector(X,1,*pb_);

#if CHEAPSIMPLE_ALGORITHM // SIMPLE and SIMPLEC but without solve, just AMG
  CheapSimple(*vx_,*px_,*vb_,*pb_);
#else

#if SIMPLER_ALGORITHM
  Simpler(*vx_,*px_,*vb_,*pb_);
#else
  Simple(*vx_,*px_,*vb_,*pb_);
#endif

#endif

  // insert solution for velocity and pressure
  mmex_.InsertVector(*vx_,0,Y);
  mmex_.InsertVector(*px_,1,Y);

  return 0;
}

/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 | taken from:                                                          |
 | Pernice, M., Tocci, M.D.:                                            |
 | A Multigrid Preconditioned Newton-Krylov method for the incomp.      |
 | Navier-Stokes equations, Siam, J. Sci. Comp. 23, pp. 398-418 (2001)  |
 *----------------------------------------------------------------------*/
void LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::Simpler(LINALG::ANA::Vector& vx, LINALG::ANA::Vector& px,
                                       LINALG::ANA::Vector& vb, LINALG::ANA::Vector& pb) const
{
  using namespace LINALG::ANA;
  SparseMatrix& A00      = (*A_)(0,0);
  SparseMatrix& A10      = (*A_)(1,0);
  SparseMatrix& A01      = (*A_)(0,1);
  SparseMatrix& diagAinv = *diagAinv_;
  SparseMatrix& S        = *S_;

  //-------------------------------------------------- L-solve / U-solve

  px = inverse(S,*psolver_,false) * ( pb - A10 * (diagAinv * vb) );

  vx = inverse(A00,*vsolver_,false) * ( vb - A01 * px );

  //------------------------------------------------ Implicit projection

  if (alpha_ != 1.0) px *= alpha_;

  *vwork2_ = diagAinv * A01 * ( inverse(S,*psolver_,false) * ( A10 * vx ) );

  vx -= vwork2_;

  return;
}

/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 | taken from:                                                          |
 | Elman, H., Howle, V.E., Shadid, J., Shuttleworth, R., Tuminaro, R.:  |
 | A taxonomy and comparison of parallel block multi-level              |
 | preconditioners for the incomp. Navier-Stokes equations.             |
 | Sandia technical report SAND2007-2761, 2007                          |
 | Also appeared in JCP                                                 |
 *----------------------------------------------------------------------*/
void LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::Simple(LINALG::ANA::Vector& vx, LINALG::ANA::Vector& px,
                                      LINALG::ANA::Vector& vb, LINALG::ANA::Vector& pb) const
{
  using namespace LINALG::ANA;
  LINALG::SparseMatrix& A00      = (*A_)(0,0);
  SparseMatrix& A10      = (*A_)(1,0);
  SparseMatrix& A01      = (*A_)(0,1);
  SparseMatrix& diagAinv = *diagAinv_;
  SparseMatrix& S        = *S_;


  //------------------------------------------------------------ L-solve

  *vwork1_ = inverse(A00,*vsolver_,false) * vb;

  px = inverse(S,*psolver_,false) * (pb - A10 * vwork1_);

  //------------------------------------------------------------ U-solve

  if (alpha_ != 1.0) px *= alpha_;

  vx = vwork1_ - diagAinv * (A01 * px);

  return;
}

/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 | is a cheaper variation from:                                         |
 | Elman, H., Howle, V.E., Shadid, J., Shuttleworth, R., Tuminaro, R.:  |
 | A taxonomy and comparison of parallel block multi-level              |
 | preconditioners for the incomp. Navier-Stokes equations.             |
 | Sandia technical report SAND2007-2761, 2007                          |
 | Also appeared in JCP                                                 |
 |                                                                      |
 | all solves replaced by single AMG sweeps                             |
 *----------------------------------------------------------------------*/
void LINALG::SOLVER::CheapSIMPLE_BlockPreconditioner::CheapSimple(LINALG::ANA::Vector& vx, LINALG::ANA::Vector& px,
                                           LINALG::ANA::Vector& vb, LINALG::ANA::Vector& pb) const
{
  SparseMatrix& A10      = (*A_)(1,0);
  SparseMatrix& A01      = (*A_)(0,1);
  SparseMatrix& diagAinv = *diagAinv_;

  //------------------------------------------------------------ L-solve
  if (vdw_)
  {
    vdwind_->Permute(&vb,&*vdwin_);
    Ppredict_->ApplyInverse(*vdwin_,*vdwout_);
    vdwind_->InvPermute(&*vdwout_,&*vwork1_);
  }
  else
    Ppredict_->ApplyInverse(vb,*vwork1_);

  *pwork1_ = pb - A10 * vwork1_;

  if (pdw_)
  {
    pdwind_->Permute(&*pwork1_,&*pdwin_);
    Pschur_->ApplyInverse(*pdwin_,*pdwout_);
    pdwind_->InvPermute(&*pdwout_,&px);
  }
  else
    Pschur_->ApplyInverse(*pwork1_,px);

  //------------------------------------------------------------ U-solve

  if (alpha_ != 1.0) px *= alpha_;

  vx = vwork1_ - diagAinv * (A01 * px);

  return;
}








