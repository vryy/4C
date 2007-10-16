/*!----------------------------------------------------------------------
\file linalg_solver.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "linalg_solver.H"

extern "C" 
{
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 02/07|
 *----------------------------------------------------------------------*/
LINALG::Solver::Solver(RefCountPtr<ParameterList> params, 
                       const Epetra_Comm& comm, FILE* outfile) :
comm_(comm),
params_(params),
outfile_(outfile),
factored_(false),
ncall_(0)
{
  // create an empty linear problem
  lp_ = rcp(new Epetra_LinearProblem());

#ifdef PARALLEL
#ifdef SPOOLES_PACKAGE
  frontmtx_      =NULL;
  newA_          =NULL;
  newY_          =NULL;
  frontETree_    =NULL;
  mtxmanager_    =NULL;
  newToOldIV_    =NULL;
  oldToNewIV_    =NULL;
  ownersIV_      =NULL;
  vtxmapIV_      =NULL;
  ownedColumnsIV_=NULL;
  solvemap_      =NULL;
  symbfacIVL_    =NULL;  
  graph_         =NULL;
  mtxY_          =NULL;
  mtxX_          =NULL;
  mtxA_          =NULL;
#endif
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 02/07|
 *----------------------------------------------------------------------*/
LINALG::Solver::~Solver()
{
#ifdef PARALLEL
#ifdef SPOOLES_PACKAGE
  if (frontmtx_)       FrontMtx_free(frontmtx_);        frontmtx_      =NULL;
  if (newA_)           InpMtx_free(newA_);              newA_          =NULL;
  if (newY_)           DenseMtx_free(newY_);            newY_          =NULL;
  if (frontETree_)     ETree_free(frontETree_);         frontETree_    =NULL;
  if (mtxmanager_)     SubMtxManager_free(mtxmanager_); mtxmanager_    =NULL;
  if (newToOldIV_)     IV_free(newToOldIV_);            newToOldIV_    =NULL;
  if (oldToNewIV_)     IV_free(oldToNewIV_);            oldToNewIV_    =NULL;
  if (ownersIV_)       IV_free(ownersIV_);              ownersIV_      =NULL;
  if (vtxmapIV_)       IV_free(vtxmapIV_);              vtxmapIV_      =NULL;
  if (ownedColumnsIV_) IV_free(ownedColumnsIV_);        ownedColumnsIV_=NULL;
  if (solvemap_)       SolveMap_free(solvemap_);        solvemap_      =NULL;
  if (graph_)          Graph_free(graph_);              graph_         =NULL;
  if (mtxY_)           DenseMtx_free(mtxY_);            mtxY_          =NULL;
  if (mtxX_)           DenseMtx_free(mtxX_);            mtxX_          =NULL;
  if (mtxA_)           InpMtx_free(mtxA_);              mtxA_          =NULL;
  if (symbfacIVL_)     IVL_free(symbfacIVL_);           symbfacIVL_    =NULL;
#endif
#endif
  return;
}

/*----------------------------------------------------------------------*
 |  reset solver (public)                                    mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Reset()
{
  A_        = null;
  Aplus_    = null;
  P_        = null;
  Pmatrix_  = null;
  x_        = null;
  b_        = null;
  lp_       = rcp(new Epetra_LinearProblem());
  factored_ = false;
  ncall_    = 0;
  amesos_   = null;
  aztec_    = null;
#ifdef PARALLEL
#ifdef SPOOLES_PACKAGE
  if (frontmtx_)       FrontMtx_free(frontmtx_);        frontmtx_      =NULL;
  if (newA_)           InpMtx_free(newA_);              newA_          =NULL;
  if (newY_)           DenseMtx_free(newY_);            newY_          =NULL;
  if (frontETree_)     ETree_free(frontETree_);         frontETree_    =NULL;
  if (mtxmanager_)     SubMtxManager_free(mtxmanager_); mtxmanager_    =NULL;
  if (newToOldIV_)     IV_free(newToOldIV_);            newToOldIV_    =NULL;
  if (oldToNewIV_)     IV_free(oldToNewIV_);            oldToNewIV_    =NULL;
  if (ownersIV_)       IV_free(ownersIV_);              ownersIV_      =NULL;
  if (vtxmapIV_)       IV_free(vtxmapIV_);              vtxmapIV_      =NULL;
  if (ownedColumnsIV_) IV_free(ownedColumnsIV_);        ownedColumnsIV_=NULL;
  if (solvemap_)       SolveMap_free(solvemap_);        solvemap_      =NULL;
  if (graph_)          Graph_free(graph_);              graph_         =NULL;
  if (mtxY_)           DenseMtx_free(mtxY_);            mtxY_          =NULL;
  if (mtxX_)           DenseMtx_free(mtxX_);            mtxX_          =NULL;
  if (mtxA_)           InpMtx_free(mtxA_);              mtxA_          =NULL;
  if (symbfacIVL_)     IVL_free(symbfacIVL_);           symbfacIVL_    =NULL;
#endif
#endif
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 02/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const LINALG::Solver& solver)
{
  solver.Print(os); 
  return os;
}

/*----------------------------------------------------------------------*
 |  print solver (public)                                    mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Print(ostream& os) const
{
  if (Comm().MyPID()==0) 
  {
    os << "============================LINALG::Solver Parameter List\n";
    os << *params_;
    os << "========================end LINALG::Solver Parameter List\n";
  }
  return;
}

/*----------------------------------------------------------------------*
 |  solve (public)                                           mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve(RefCountPtr<Epetra_CrsMatrix> matrix,
                           RefCountPtr<Epetra_Vector>    x,
                           RefCountPtr<Epetra_Vector>    b,
                           bool refactor,
                           bool reset)
{
  // reset data flags
  if (reset) 
  {
    Reset();
    refactor = true;
  }
  
  // set the data passed to the method
  if (refactor) A_ = matrix;
  x_ = x;
  b_ = b;
  
  // set flag indicating that problem should be refactorized
  if (refactor) factored_ = false;
  
  // fill the linear problem
  lp_->SetRHS(b_.get());
  lp_->SetLHS(x_.get());
  lp_->SetOperator(A_.get());
  
  // decide what solver to use
  string solvertype = Params().get("solver","none");
  if      ("lapack" ==solvertype) 
    Solve_lapack(reset);
  else if ("klu"    ==solvertype) 
    Solve_klu(reset);
  else if ("umfpack"==solvertype) 
    Solve_umfpack(reset);
#ifdef PARALLEL
  else if ("superlu"==solvertype) 
    Solve_superlu(reset);
#endif
#ifdef PARALLEL
#ifdef SPOOLES_PACKAGE
  else if ("spooles"==solvertype) 
    Solve_spooles(reset);
#endif
#endif
  else if ("aztec"  ==solvertype) 
    Solve_aztec(reset);
  else if ("none"   ==solvertype) 
    dserror("Unknown type of solver");
  
  factored_ = true;
  ncall_++;
  
  return;
}

/*----------------------------------------------------------------------*
 |  solve (public)                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve(RefCountPtr<Epetra_Operator>  Operator,
                           RefCountPtr<Epetra_Vector>    x,
                           RefCountPtr<Epetra_Vector>    b,
                           bool refactor,
                           bool reset)
{
  // reset data flags
  if (reset) 
  {
    Reset();
    refactor = true;
  }
  
  // set the data passed to the method
  if (refactor) A_ = Operator;
  x_ = x;
  b_ = b;
  
  // set flag indicating that problem should be refactorized
  if (refactor) factored_ = false;
  
  // fill the linear problem
  lp_->SetRHS(b_.get());
  lp_->SetLHS(x_.get());
  lp_->SetOperator(A_.get());
  
  // decide what solver to use
  string solvertype = Params().get("solver","none");
  if      ("aztec" ==solvertype) 
    Solve_aztec(reset);
  else if ("none"   ==solvertype) 
    dserror("Unknown type of solver");
  else
    dserror("Unsupported type of solver");
  
  
  factored_ = true;
  ncall_++;
  
  return;
}

/*----------------------------------------------------------------------*
 |  solve vm3 (public)                                     vgravem 06/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve(RefCountPtr<Epetra_CrsMatrix> matrix,
                           RefCountPtr<Epetra_CrsMatrix> matrix2,
                           RefCountPtr<Epetra_Vector>    x,
                           RefCountPtr<Epetra_Vector>    b,
                           bool refactor,
                           bool reset)
{
  // see whether Operator is a Epetra_CrsMatrix
  Epetra_CrsMatrix* tmp = dynamic_cast<Epetra_CrsMatrix*>(A_.get());
  if (!tmp) dserror("vm3 only with Epetra_Operator being an Epetra_CrsMatrix!");
  RCP<Epetra_CrsMatrix> A = rcp(tmp);
  A.release();
  
  // reset data flags
  if (reset) 
  {
    Reset();
    refactor = true;
  }

  
  // set the data passed to the method
  if (refactor)
  {
    Aplus_ = matrix;
    A_     = matrix2;
  }
  x_ = x;
  b_ = b;
  
  // set flag indicating that problem should be refactorized
  if (refactor) factored_ = false;
  
  // fill the linear problem
  lp_->SetRHS(b_.get());
  lp_->SetLHS(x_.get());
  lp_->SetOperator(Aplus_.get());
  lp_->SetOperator(A_.get());
  
  // decide what solver to use
  string solvertype = Params().get("solver","none");
  // if vm3 was not selected, return to usual solver
  if ("vm3"   !=solvertype)
  {
    LINALG::Solver::Solve(matrix,x,b,refactor,reset);
    return;
  }

  // extract the ML parameters and initialize the solver
  ParameterList&  mllist = Params().sublist("ML Parameters");
  // cout << "Parameter list:\n" << mllist;
  vm3_solver_ = rcp(new VM3_Solver::VM3_Solver(Aplus_, A, mllist,true) );

  // Apply the solver. Convergence check and iteration have to be
  // performed from outside, since Solve is a nonlinear FAS scheme
  B_ = b_.get();
  X_ = x_.get();
  int err = vm3_solver_->VM3_Solver::Solve(*B_,*X_);
  if (err) dserror("VM3_Solver::Solve returned an err");

  ncall_++;
  
  return;
}


/*----------------------------------------------------------------------*
 |  solve (protected)                                        mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve_aztec(const bool reset)
{
  if (!Params().isSublist("Aztec Parameters")) 
    dserror("Do not have aztec parameter list");
  ParameterList& azlist = Params().sublist("Aztec Parameters");

  // see whether Operator is a Epetra_CrsMatrix
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(A_.get());

  // decide whether we recreate preconditioners
  bool create = false;
  int  reuse  = azlist.get("reuse",0);
  if      (reset)            create = true;
  else if (!Ncall())         create = true;
  else if (!reuse)           create = true;
  else if (Ncall()%reuse==0) create = true;

  // Allocate an aztec solver
  if (create)
  {
    // create an aztec solver
    aztec_ = rcp(new AztecOO());
    aztec_->SetAztecDefaults();
    //aztec_->SetParameters(azlist,false); moved this further down gee
  }

  // decide whether we do what kind of scaling
  bool scaling_infnorm = false;
  bool scaling_symdiag = false;
  string scaling = azlist.get("scaling","none");
  if (scaling=="none");
  else if (scaling=="infnorm")
  {
    scaling_infnorm = true;
    scaling_symdiag = false;
  }
  else if (scaling=="symmetric")
  {
    scaling_infnorm = false;
    scaling_symdiag = true;
  }
  else dserror("Unknown type of scaling found in parameter list");	

  if (!A)
  {
    scaling_infnorm = false;
    scaling_symdiag = false;
  }

  // do infnorm scaling
  RefCountPtr<Epetra_Vector> rowsum;
  RefCountPtr<Epetra_Vector> colsum;
  if (scaling_infnorm)
  {
    rowsum = rcp(new Epetra_Vector(A->RowMap(),false));
    colsum = rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*rowsum);
    A->InvColSums(*colsum);
    lp_->LeftScale(*rowsum);
    lp_->RightScale(*colsum);
  }
  
  // do symmetric diagonal scaling
  RefCountPtr<Epetra_Vector> diag;
  if (scaling_symdiag)
  {
    Epetra_Vector invdiag(A->RowMap(),false);
    diag = rcp(new Epetra_Vector(A->RowMap(),false));
    A->ExtractDiagonalCopy(*diag);
    invdiag.Reciprocal(*diag);
    lp_->LeftScale(invdiag);
    lp_->RightScale(invdiag);
  }
  
  // pass linear problem to aztec
  aztec_->SetProblem(*lp_);
  // don't want linear problem to alter our atzec parameters (idiot feature!)
  aztec_->SetParameters(azlist,false);

  // get type of preconditioner and build either Ifpack or ML
  // if we have an ifpack parameter list, we do ifpack
  // if we have an ml parameter list we do ml
  bool doifpack = Params().isSublist("IFPACK Parameters");
  bool doml     = Params().isSublist("ML Parameters");
  if (!A)
  {
    doifpack = false;
    doml     = false;
  }


  // do ifpack if desired
  if (create && doifpack)
  {
    ParameterList&  ifpacklist = Params().sublist("IFPACK Parameters");
    // create a copy of the scaled matrix
    // so we can reuse the precondition
    Pmatrix_ = rcp(new Epetra_CrsMatrix(*A));
    // get the type of ifpack preconditioner from aztec
    string prectype = azlist.get("preconditioner","ILU");
    int    overlap  = azlist.get("AZ_overlap",0);
    Ifpack Factory;
    Ifpack_Preconditioner* prec = Factory.Create(prectype,Pmatrix_.get(),overlap);
    prec->SetParameters(ifpacklist);
    prec->Initialize();
    prec->Compute();
    P_ = rcp(prec);
  }
  
  // do ml if desired
  if (create && doml)
  {
    ParameterList&  mllist = Params().sublist("ML Parameters");
    // create a copy of the scaled matrix
    // so we can reuse the precondition
    Pmatrix_ = rcp(new Epetra_CrsMatrix(*A));
    P_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*Pmatrix_,mllist,true));
  }
  
  if (doifpack || doml)
    aztec_->SetPrecOperator(P_.get());

  // iterate on the solution
  int iter = azlist.get("AZ_max_iter",500);
  double tol = azlist.get("AZ_tol",1.0e-6);
  aztec_->Iterate(iter,tol);
  
  // check status of solution process
  const double* status = aztec_->GetAztecStatus();
  if (status[AZ_why] != AZ_normal)
  {
    bool resolve = false;
    if (status[AZ_why] == AZ_breakdown)
    {
      if (Comm().MyPID()==0)
        printf("Numerical breakdown in AztecOO, try again with KLU/SuperLU\n");
      resolve = true;
    }
    else if (status[AZ_why] == AZ_ill_cond)
    {
      if (Comm().MyPID()==0)
        printf("Problem is near singular in AztecOO, try again with KLU/SuperLU\n");
      resolve = true;
    }
    else if (status[AZ_why] == AZ_maxits)
    {
      if (Comm().MyPID()==0)
        printf("Max iterations reached in AztecOO, try again with KLU/SuperLU\n");
      resolve = true;
    }
    if (resolve)
    {
#if 0
#ifdef PARALLEL
      if (Comm().MyPID()==0) cout << "Retrying using SuperLU\n"; fflush(stdout);
      Amesos_Superludist superlusolver(*lp_);
      int err = superlusolver.SymbolicFactorization();
      if (err) dserror("SuperLU.SymbolicFactorization() returned %d",err);
      err     = superlusolver.NumericFactorization();
      if (err) dserror("SuperLU.NumericFactorization() returned %d",err);
      err     = superlusolver.Solve();
      if (err) dserror("SuperLU.Solve() returned %d",err);
#else
      if (Comm().MyPID()==0) cout << "Retrying using KLU\n"; fflush(stdout);
      Amesos_Klu klusolver(*lp_);
      int err = klusolver.SymbolicFactorization();
      if (err) dserror("Amesos_Klu.SymbolicFactorization() returned %d",err);
      err     = klusolver.NumericFactorization();
      if (err) dserror("Amesos_Klu.NumericFactorization() returned %d",err);
      err     = klusolver.Solve();
      if (err) dserror("Amesos_Klu.Solve() returned %d",err);        
#endif
#endif      
    }
  } // if (status[AZ_why] != AZ_normal)

  // undo scaling
  if (scaling_infnorm)
  {
    Epetra_Vector invrowsum(A->RowMap(),false);
    invrowsum.Reciprocal(*rowsum);
    rowsum = null;
    Epetra_Vector invcolsum(A->RowMap(),false);
    invcolsum.Reciprocal(*colsum);
    colsum = null;
    lp_->LeftScale(invrowsum);
    lp_->RightScale(invcolsum);
  }
  if (scaling_symdiag)
  {
    lp_->LeftScale(*diag);
    lp_->RightScale(*diag);
    diag = null;
  }
  


  // print some output if desired
  if (Comm().MyPID()==0 && outfile_)
  {
    fprintf(outfile_,"AztecOO: unknowns/iterations/time %d  %d  %f\n",
            A_->OperatorRangeMap().NumGlobalElements(),(int)status[AZ_its],status[AZ_solve_time]);
    fflush(outfile_);
  }

  return;
}



/*----------------------------------------------------------------------*
 |  solve (protected)                                        mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve_superlu(const bool reset)
{
#ifdef PARALLEL
  if (reset || !IsFactored())
    amesos_ = rcp(new Amesos_Superludist(*lp_));

  if (amesos_==null) dserror("No solver allocated");

  // Problem has not been factorized before
  if (!IsFactored())
  {
    int err = amesos_->SymbolicFactorization();
    if (err) dserror("Amesos::SymbolicFactorization returned an err");
    err = amesos_->NumericFactorization();
    if (err) dserror("Amesos::NumericFactorization returned an err");
  }
  
  int err = amesos_->Solve();
  if (err) dserror("Amesos::Solve returned an err");
#else
  dserror("Distributed SuperLU only in parallel");
#endif    //! system of equations
  RefCountPtr<Epetra_CrsMatrix>     A_;       

  return;
}

/*----------------------------------------------------------------------*
 |  solve (protected)                                        mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve_umfpack(const bool reset)
{
  if (reset || !IsFactored())
    amesos_ = rcp(new Amesos_Umfpack(*lp_));

  if (amesos_==null) dserror("No solver allocated");

  // Problem has not been factorized before
  if (!IsFactored())
  {
    bool symmetric = Params().get("symmetric",false);
    amesos_->SetUseTranspose(symmetric);
    int err = amesos_->SymbolicFactorization();
    if (err) dserror("Amesos::SymbolicFactorization returned an err");
    err = amesos_->NumericFactorization();
    if (err) dserror("Amesos::NumericFactorization returned an err");
  }
  
  int err = amesos_->Solve();
  if (err) dserror("Amesos::Solve returned an err");
  
  return;
}

/*----------------------------------------------------------------------*
 |  solve (protected)                                        mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve_klu(const bool reset)
{
  if (reset || !IsFactored())
    amesos_ = rcp(new Amesos_Klu(*lp_));

  if (amesos_==null) dserror("No solver allocated");

  // Problem has not been factorized before
  if (!IsFactored())
  {
    bool symmetric = Params().get("symmetric",false);
    amesos_->SetUseTranspose(symmetric);
    int err = amesos_->SymbolicFactorization();
    if (err) dserror("Amesos::SymbolicFactorization returned an err");
    err = amesos_->NumericFactorization();
    if (err) dserror("Amesos::NumericFactorization returned an err");
  }
  
  int err = amesos_->Solve();
  if (err) dserror("Amesos::Solve returned an err");
  
  return;
}

/*----------------------------------------------------------------------*
 |  solve (protected)                                        mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve_lapack(const bool reset)
{
  if (reset || !IsFactored())
    amesos_ = rcp(new Amesos_Lapack(*lp_));

  if (amesos_==null) dserror("No solver allocated");

  // Problem has not been factorized before
  if (!IsFactored())
  {
    int err = amesos_->SymbolicFactorization();
    if (err) dserror("Amesos::SymbolicFactorization returned an err");
    err = amesos_->NumericFactorization();
    if (err) dserror("Amesos::NumericFactorization returned an err");
  }
  
  int err = amesos_->Solve();
  if (err) dserror("Amesos::Solve returned an err");
  
  return;
}

/*----------------------------------------------------------------------*
 |  translate solver parameters (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::TranslateSolverParameters(ParameterList& params, 
                                          struct _SOLVAR* actsolv) const
{
  // switch type of solver
  switch (actsolv->solvertyp)
  {
#ifdef PARALLEL
  case superlu://============================== superlu solver (parallel only)
    params.set("solver","superlu");
    params.set("symmetric",false);
  break;
#endif
  case amesos_klu_sym://====================================== Tim Davis' KLU
    params.set("solver","klu");
    params.set("symmetric",true);
  break;
  case amesos_klu_nonsym://=================================== Tim Davis' KLU
    params.set("solver","klu");
    params.set("symmetric",false);
  break;
  case umfpack://========================================= Tim Davis' Umfpack
    params.set("solver","umfpack");
    params.set("symmetric",false);
  break;
  case lapack_sym://================================================== Lapack
    params.set("solver","lapack");
    params.set("symmetric",true);
  break;
  case lapack_nonsym://=============================================== Lapack
    params.set("solver","lapack");
    params.set("symmetric",false);
  break;
  case vm3://=========================================================== VM3
  {
    params.set("solver","vm3");
    ParameterList& mllist = params.sublist("ML Parameters");
    ML_Epetra::SetDefaults("SA",mllist);
    AZVAR* azvar = actsolv->azvar;
    mllist.set("output"                          ,azvar->mlprint);
    if (azvar->mlprint==10)
      mllist.set("print unused"                  ,1);
    else
    mllist.set("print unused"                  ,-2);
    mllist.set("increasing or decreasing"        ,"increasing");
    mllist.set("coarse: max size"                ,azvar->mlcsize);
    mllist.set("max levels"                      ,azvar->mlmaxlevel);
    mllist.set("smoother: pre or post"           ,"both");
    mllist.set("aggregation: use tentative restriction",true);
    mllist.set("aggregation: threshold"          ,azvar->ml_threshold);
    mllist.set("aggregation: damping factor"     ,azvar->mldamp_prolong);
    mllist.set("aggregation: nodes per aggregate",azvar->mlaggsize);
    switch (azvar->mlcoarsentype)
    {
      case 0:  mllist.set("aggregation: type","Uncoupled");  break;
      case 1:  mllist.set("aggregation: type","METIS");      break;
      case 2:  mllist.set("aggregation: type","VBMETIS");    break;
      case 3:  mllist.set("aggregation: type","MIS");        break;
      default: dserror("Unknown type of coarsening for ML"); break;
    }
    // set ml smoothers -> write ml list
    for (int i=0; i<azvar->mlmaxlevel-1; ++i)
    {
      char levelstr[11];
      sprintf(levelstr,"(level %d)",i);
      int type;
      double damp;
      if (i==0)
      {
        type = azvar->mlsmotype_fine;
        damp = azvar->mldamp_fine;
      }
      else if (i < azvar->mlmaxlevel-1)
      {
        type = azvar->mlsmotype_med;
        damp = azvar->mldamp_med;
      }
      else
      {
        type = azvar->mlsmotype_coarse;
        damp = azvar->mldamp_coarse;
      }
      switch (type)
      {
      case 0:
        mllist.set("smoother: type"                   ,"symmetric Gauss-Seidel");
        mllist.set("smoother: sweeps"                 ,azvar->mlsmotimes[i]);
        mllist.set("smoother: damping factor"         ,damp);
      break;
      case 1:
        mllist.set("smoother: type"                    ,"Jacobi");
        mllist.set("smoother: sweeps"                  ,azvar->mlsmotimes[i]);
        mllist.set("smoother: damping factor"          ,damp);
      break;
      case 2:
        mllist.set("smoother: type "+(string)levelstr                    ,"MLS");
        mllist.set("smoother: MLS polynomial order "+(string)levelstr    ,azvar->mlsmotimes[i]);
      break;
      case 3:
        mllist.set("smoother: type (level 0)"                            ,"MLS");
        mllist.set("smoother: MLS polynomial order "+(string)levelstr    ,-azvar->mlsmotimes[i]);
      break;
      case 4:
        mllist.set("smoother: type "+(string)levelstr                    ,"IFPACK");
        mllist.set("smoother: ifpack type "+(string)levelstr             ,"ILU");
        mllist.set("smoother: ifpack overlap "+(string)levelstr          ,0);
        mllist.sublist("smoother: ifpack list").set("fact: level-of-fill",azvar->mlsmotimes[i]);
        mllist.sublist("smoother: ifpack list").set("schwarz: reordering type","rcm");
      break;
      case 5:
        mllist.set("smoother: type "+(string)levelstr,"Amesos-KLU");
      break;
#ifdef PARALLEL
      case 6:
        mllist.set("smoother: type "+(string)levelstr,"Amesos-Superludist");
      break;
#endif
      default: dserror("Unknown type of smoother for ML"); break;
      } // switch (type)
    } // for (int i=0; i<azvar->mlmaxlevel-1; ++i)
    // set coarse grid solver
    const int coarse = azvar->mlmaxlevel-1;
    switch (azvar->mlsmotype_coarse)
    {
      case 0:
        mllist.set("coarse: type"          ,"symmetric Gauss-Seidel");
        mllist.set("coarse: sweeps"        , azvar->mlsmotimes[coarse]);
        mllist.set("coarse: damping factor",azvar->mldamp_coarse);
      break;
      case 1:
        mllist.set("coarse: type"          ,"Jacobi");
        mllist.set("coarse: sweeps"        , azvar->mlsmotimes[coarse]);
        mllist.set("coarse: damping factor",azvar->mldamp_coarse);
      break;
      case 2:
        mllist.set("coarse: type"                ,"MLS");
        mllist.set("coarse: MLS polynomial order",azvar->mlsmotimes[coarse]);
      break;
      case 3:
        mllist.set("coarse: type"                ,"MLS");
        mllist.set("coarse: MLS polynomial order",-azvar->mlsmotimes[coarse]);
      break;
      case 4:
        mllist.set("coarse: type"          ,"IFPACK");
        mllist.set("coarse: ifpack type"   ,"ILU");
        mllist.set("coarse: ifpack overlap",0);
        mllist.sublist("coarse: ifpack list").set("fact: level-of-fill",azvar->mlsmotimes[coarse]);
        mllist.sublist("coarse: ifpack list").set("schwarz: reordering type","rcm");
      break;
      case 5:
        mllist.set("coarse: type","Amesos-KLU");
      break;
      case 6:
        mllist.set("coarse: type","Amesos-Superludist");
      break;
      default: dserror("Unknown type of coarse solver for ML"); break;
    } // switch (azvar->mlsmotype_coarse)
    // default values for nullspace
    mllist.set("PDE equations",1);
    mllist.set("null space: dimension",1);
    mllist.set("null space: type","pre-computed");
    mllist.set("null space: add default vectors",false);
    mllist.set<double*>("null space: vectors",NULL);
  }
  break;
  case aztec_msr://================================================= AztecOO
  {
    params.set("solver","aztec");
    params.set("symmetric",false);
    AZVAR* azvar = actsolv->azvar;
    ParameterList& azlist = params.sublist("Aztec Parameters");
    //--------------------------------- set scaling of linear problem
    if (azvar->azscal==1)
      azlist.set("scaling","symmetric");
    else if (azvar->azscal==2)
      azlist.set("scaling","infnorm");
    else
      azlist.set("scaling","none");
    //--------------------------------------------- set type of solver
    switch (azvar->azsolvertyp)
    {
    case azsolv_CG:       azlist.set("AZ_solver",AZ_cg);       break;
    case azsolv_GMRES:    azlist.set("AZ_solver",AZ_gmres);    break;
    case azsolv_CGS:      azlist.set("AZ_solver",AZ_cgs);      break;
    case azsolv_BiCGSTAB: azlist.set("AZ_solver",AZ_bicgstab); break;
    case azsolv_LU:       azlist.set("AZ_solver",AZ_lu);       break;
    case azsolv_TFQMR:    azlist.set("AZ_solver",AZ_tfqmr);    break;
    default: dserror("Unknown solver for AztecOO");            break;
    }
    //------------------------------------- set type of preconditioner
    switch (azvar->azprectyp)
    {
    case azprec_none:
      azlist.set("AZ_precond",AZ_none);
      azlist.set("AZ_subdomain_solve",AZ_none);
      azlist.set("preconditioner",AZ_none);
    break;
    case azprec_ILUT:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("preconditioner","ILUT");
    break;
    case azprec_ILU:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("preconditioner","ILU");
    break;
    case azprec_Jacobi:
      azlist.set("AZ_precond",AZ_Jacobi);
    break;
    case azprec_Neumann:
      azlist.set("AZ_precond",AZ_Neumann);
    break;
    case azprec_Least_Squares:
      azlist.set("AZ_precond",AZ_ls);
    break;
    case azprec_SymmGaussSeidel:
      azlist.set("AZ_precond",AZ_sym_GS);
    break;
    case azprec_LU:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("preconditioner","Amesos");
    break;
    case azprec_RILU:
      azlist.set("AZ_precond",AZ_dom_decomp);
      azlist.set("AZ_subdomain_solve",AZ_rilu);
      azlist.set("AZ_graph_fill",azvar->azgfill);
    break;
    case azprec_ICC:
      // using ifpack
      azlist.set("AZ_precond",AZ_user_precond);
      azlist.set("preconditioner","IC");
    break;
    case azprec_ML:
    case azprec_MLfluid:
    case azprec_MLfluid2:
      azlist.set("AZ_precond",AZ_user_precond);
    break;
    default:
      dserror("Unknown preconditioner for AztecOO");
    break;
    }
    //------------------------------------- set other aztec parameters
    azlist.set("AZ_kspace",azvar->azsub);
    azlist.set("AZ_max_iter",azvar->aziter);
    azlist.set("AZ_overlap",0);
    azlist.set("AZ_type_overlap",AZ_symmetric);
    azlist.set("AZ_poly_ord",azvar->azpoly);
    if (!azvar->azoutput)
      azlist.set("AZ_output",AZ_none);             // AZ_none AZ_all AZ_warnings AZ_last 10
    else
      azlist.set("AZ_output",azvar->azoutput);
    azlist.set("AZ_diagnostics",AZ_none);          // AZ_none AZ_all
    azlist.set("AZ_conv",azvar->azconv);
    azlist.set("AZ_tol",azvar->aztol);
    azlist.set("AZ_drop",azvar->azdrop);
    azlist.set("AZ_scaling",AZ_none);              
    azlist.set("AZ_keep_info",0);
    // set reuse parameters
    azlist.set("ncall",0);                         // counting number of solver calls
    azlist.set("reuse",azvar->azreuse);            // reuse info for n solver calls
    //-------------------------------- set parameters for Ifpack if used
    if (azvar->azprectyp == azprec_ILU  ||
        azvar->azprectyp == azprec_ILUT ||
        azvar->azprectyp == azprec_ICC  ||
        azvar->azprectyp == azprec_LU   )
    {
      ParameterList& ifpacklist = params.sublist("IFPACK Parameters");
      ifpacklist.set("fact: drop tolerance",azvar->azdrop);
      ifpacklist.set("fact: level-of-fill",azvar->azgfill);
      ifpacklist.set("fact: ilut level-of-fill",azvar->azfill);
      ifpacklist.set("schwarz: combine mode","Add"); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
      ifpacklist.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
      ifpacklist.set("amesos: solver type", "Amesos_Klu"); // can be "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"
    }
    //------------------------------------- set parameters for ML if used
    if (azvar->azprectyp == azprec_ML      ||
        azvar->azprectyp == azprec_MLfluid ||
        azvar->azprectyp == azprec_MLfluid2 )
    {
      ParameterList& mllist = params.sublist("ML Parameters");
      ML_Epetra::SetDefaults("SA",mllist);
      switch (azvar->azprectyp)
      {
      case azprec_ML: // do nothing, this is standard
      break;
      case azprec_MLfluid: // unsymmetric, unsmoothed restruction
        mllist.set("aggregation: use tentative restriction",true);
      break;
      case azprec_MLfluid2: // full Pretrov-Galerkin unsymmetric smoothed
        mllist.set("energy minimization: enable",true);
        mllist.set("energy minimization: type",3); // 1,2,3 cheap -> expensive
        mllist.set("aggregation: block scaling",false); 
      break;
      default: dserror("Unknown type of ml preconditioner");
      }
      mllist.set("output"                          ,azvar->mlprint);
      if (azvar->mlprint==10)
        mllist.set("print unused"                  ,1);
      else
        mllist.set("print unused"                  ,-2);
      mllist.set("increasing or decreasing"        ,"increasing");
      mllist.set("coarse: max size"                ,azvar->mlcsize);
      mllist.set("max levels"                      ,azvar->mlmaxlevel);
      mllist.set("smoother: pre or post"           ,"both");
      mllist.set("aggregation: threshold"          ,azvar->ml_threshold);
      mllist.set("aggregation: damping factor"     ,azvar->mldamp_prolong);
      mllist.set("aggregation: nodes per aggregate",azvar->mlaggsize);
      switch (azvar->mlcoarsentype)
      {
        case 0:  mllist.set("aggregation: type","Uncoupled");  break;
        case 1:  mllist.set("aggregation: type","METIS");      break;
        case 2:  mllist.set("aggregation: type","VBMETIS");    break;
        case 3:  mllist.set("aggregation: type","MIS");        break;
        default: dserror("Unknown type of coarsening for ML"); break;
      }
      // set ml smoothers
      for (int i=0; i<azvar->mlmaxlevel-1; ++i)
      {
        char levelstr[11];
        sprintf(levelstr,"(level %d)",i);
        int type;
        double damp;
        if (i==0)
        {
          type = azvar->mlsmotype_fine;
          damp = azvar->mldamp_fine;
        }
        else if (i < azvar->mlmaxlevel-1)
        {
          type = azvar->mlsmotype_med;
          damp = azvar->mldamp_med;
        }
        else
        {
          type = azvar->mlsmotype_coarse;
          damp = azvar->mldamp_coarse;
        }
        switch (type)
        {
        case 0:
          mllist.set("smoother: type "+(string)levelstr                    ,"symmetric Gauss-Seidel");
          mllist.set("smoother: sweeps "+(string)levelstr                  ,azvar->mlsmotimes[i]);
          mllist.set("smoother: damping factor "+(string)levelstr          ,damp);
        break;
        case 1:
          mllist.set("smoother: type "+(string)levelstr                    ,"Jacobi");
          mllist.set("smoother: sweeps "+(string)levelstr                  ,azvar->mlsmotimes[i]);
          mllist.set("smoother: damping factor "+(string)levelstr          ,damp);
        break;
        case 2:
          mllist.set("smoother: type "+(string)levelstr                    ,"MLS");
          mllist.set("smoother: MLS polynomial order "+(string)levelstr    ,azvar->mlsmotimes[i]);
        break;
        case 3:
          mllist.set("smoother: type (level 0)"                            ,"MLS");
          mllist.set("smoother: MLS polynomial order "+(string)levelstr    ,-azvar->mlsmotimes[i]);
        break;
        case 4:
          mllist.set("smoother: type "+(string)levelstr                    ,"IFPACK");
          mllist.set("smoother: ifpack type "+(string)levelstr             ,"ILU");
          mllist.set("smoother: ifpack overlap "+(string)levelstr          ,0);
          mllist.sublist("smoother: ifpack list").set("fact: level-of-fill",azvar->mlsmotimes[i]);
          mllist.sublist("smoother: ifpack list").set("schwarz: reordering type","rcm");
        break;
        case 5:
          mllist.set("smoother: type "+(string)levelstr,"Amesos-KLU");
        break;
#ifdef PARALLEL
        case 6:
          mllist.set("smoother: type "+(string)levelstr,"Amesos-Superludist");
        break;
#endif
        default: dserror("Unknown type of smoother for ML"); break;
        } // switch (type)
      } // for (int i=0; i<azvar->mlmaxlevel-1; ++i)
      // set coarse grid solver
      const int coarse = azvar->mlmaxlevel-1;
      switch (azvar->mlsmotype_coarse)
      {
        case 0:
          mllist.set("coarse: type"          ,"symmetric Gauss-Seidel");
          mllist.set("coarse: sweeps"        , azvar->mlsmotimes[coarse]);
          mllist.set("coarse: damping factor",azvar->mldamp_coarse);
        break;
        case 1:
          mllist.set("coarse: type"          ,"Jacobi");
          mllist.set("coarse: sweeps"        , azvar->mlsmotimes[coarse]);
          mllist.set("coarse: damping factor",azvar->mldamp_coarse);
        break;
        case 2:
          mllist.set("coarse: type"                ,"MLS");
          mllist.set("coarse: MLS polynomial order",azvar->mlsmotimes[coarse]);
        break;
        case 3:
          mllist.set("coarse: type"                ,"MLS");
          mllist.set("coarse: MLS polynomial order",-azvar->mlsmotimes[coarse]);
        break;
        case 4:
          mllist.set("coarse: type"          ,"IFPACK");
          mllist.set("coarse: ifpack type"   ,"ILU");
          mllist.set("coarse: ifpack overlap",0);
          mllist.sublist("coarse: ifpack list").set("fact: level-of-fill",azvar->mlsmotimes[coarse]);
          mllist.sublist("coarse: ifpack list").set("schwarz: reordering type","rcm");
        break;
        case 5:
          mllist.set("coarse: type","Amesos-KLU");
        break;
        case 6:
          mllist.set("coarse: type","Amesos-Superludist");
        break;
        default: dserror("Unknown type of coarse solver for ML"); break;
      } // switch (azvar->mlsmotype_coarse)
      // default values for nullspace
      mllist.set("PDE equations",1);
      mllist.set("null space: dimension",1);
      mllist.set("null space: type","pre-computed");
      mllist.set("null space: add default vectors",false);
      mllist.set<double*>("null space: vectors",NULL);
    } // if ml preconditioner
  }	
  break;
#ifdef PARALLEL
  case SPOOLES_sym://================================== Spooles (parallel only)
  case SPOOLES_nonsym:
    params.set("solver","spooles");
    params.set("symmetric",false);
  break;
#endif
  default:
    dserror("Unsupported type of solver");
  break;
  }
  return;
}




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
