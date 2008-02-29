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

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "linalg_solver.H"
#include "linalg_mlapi_operator.H"
#include "simpler_operator.H"
#include "linalg_systemmatrix.H"

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
 |  adapt tolerance (public)                                 mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::Solver::AdaptTolerance(const double desirednlnres,
                                    const double currentnlnres,
                                    const double better)
{
  if (!Params().isSublist("Aztec Parameters")) return;
  const int myrank = Comm().MyPID();
  ParameterList& azlist = Params().sublist("Aztec Parameters");
  int output   = azlist.get<int>("AZ_output",1);
  int convtest = azlist.get<int>("AZ_conv",AZ_noscaled);
  if (convtest != AZ_r0) dserror("Using convergence adaptivity: Use AZ_r0 in input file");
  bool havesavedvalue = azlist.isParameter("AZ_tol save");
  if (!havesavedvalue)
  {
    if (!azlist.isParameter("AZ_tol"))
    {
      cout << azlist;
      dserror("No Aztec tolerance in ParameterList");
    }
    azlist.set<double>("AZ_tol save",azlist.get<double>("AZ_tol",1.e-8));
  }
  double tol = azlist.get<double>("AZ_tol save",1.e-8);
  if (!myrank && output)
    printf("                --- Aztec input   relative tolerance %10.3E\n",tol);
  if (currentnlnres*tol < desirednlnres)
  {
    double tolnew = desirednlnres*better/currentnlnres;
    if (tolnew<tol) tolnew = tol;
    if (!myrank && output && tolnew > tol)
      printf("                *** Aztec adapted relative tolerance %10.3E\n",tolnew);
    azlist.set<double>("AZ_tol",tolnew);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  adapt tolerance (public)                                 mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::Solver::ResetTolerance()
{
  if (!Params().isSublist("Aztec Parameters")) return;
  ParameterList& azlist = Params().sublist("Aztec Parameters");
  bool havesavedvalue = azlist.isParameter("AZ_tol save");
  if (!havesavedvalue) return;
  azlist.set<double>("AZ_tol",azlist.get<double>("AZ_tol save",1.e-8));
  return;
}

/*----------------------------------------------------------------------*
 |  solve (public)                                           mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve(RefCountPtr<Epetra_Operator>  matrix,
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
  if ("aztec"  ==solvertype)
    Solve_aztec(reset);
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
  else if ("lapack" ==solvertype)
    Solve_lapack(reset);
  else if ("none"   ==solvertype)
    dserror("Unknown type of solver");

  factored_ = true;
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

  // Allocate an aztec solver with default parameters
  if (create)
  {
    // create an aztec solver
    aztec_ = rcp(new AztecOO());
    aztec_->SetAztecDefaults();
    // tell aztec to which stream to write
    aztec_->SetOutputStream(std::cout);
    aztec_->SetErrorStream(std::cerr);
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

  // Don't want linear problem to alter our aztec parameters (idiot feature!)
  // this is why we set our list here AFTER the linear problem has been set
  //
  // We don't want to use Aztec's scaling capabilities as we prefer to do
  // the scaling ourselves (so we precisely know what happens)
  // Therefore set scaling parameter to none and reset it after aztec has made
  // its internal copy of the parameter list
  azlist.set("scaling","none");
  aztec_->SetParameters(azlist,false);
  azlist.set("scaling",scaling);

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
    // see whether we use standard ml or our own mlapi operator
    const bool domlapioperator = mllist.get<bool>("LINALG::AMG_Operator",false);
    // see whether we have a SIMPLER sublist
    const bool dosimpler = Params().isSublist("SIMPLER");
    if (domlapioperator)
    {
      // create a copy of the scaled matrix
      // so we can reuse the preconditioner several times
      Pmatrix_ = rcp(new Epetra_CrsMatrix(*A));
      P_ = rcp(new LINALG::AMG_Operator(Pmatrix_,mllist,true));
    }
    else if (dosimpler)
    {
      // SIMPLER does not need copy of preconditioning matrix to live
      RCP<Epetra_CrsMatrix> rcpA = rcp(A,false);
      P_ = rcp(new LINALG::SIMPLER_Operator(rcpA,
                                            Params(),
                                            Params().sublist("SIMPLER"),
                                            outfile_));
      Pmatrix_ = null;
    }
    else
    {
      // create a copy of the scaled matrix
      // so we can reuse the preconditioner several times
      Pmatrix_ = rcp(new Epetra_CrsMatrix(*A));
      P_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*Pmatrix_,mllist,true));
      // for debugging ML
      //dynamic_cast<ML_Epetra::MultiLevelPreconditioner&>(*P_).PrintUnused(0);
    }
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
        printf("Numerical breakdown in AztecOO\n");
      resolve = true;
    }
    else if (status[AZ_why] == AZ_ill_cond)
    {
      if (Comm().MyPID()==0)
        printf("Problem is near singular in AztecOO\n");
      resolve = true;
    }
    else if (status[AZ_why] == AZ_maxits)
    {
      if (Comm().MyPID()==0)
        printf("Max iterations reached in AztecOO\n");
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
  {
    reindexer_ = rcp(new EpetraExt::LinearProblem_Reindex(NULL));
    amesos_ = rcp(new Amesos_Superludist((*reindexer_)(*lp_)));
  }

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
  {
    reindexer_ = rcp(new EpetraExt::LinearProblem_Reindex(NULL));
    amesos_ = rcp(new Amesos_Umfpack((*reindexer_)(*lp_)));
  }

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
  {
    reindexer_ = rcp(new EpetraExt::LinearProblem_Reindex(NULL));
    amesos_ = rcp(new Amesos_Klu((*reindexer_)(*lp_)));
  }

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
  {
    reindexer_ = rcp(new EpetraExt::LinearProblem_Reindex(NULL));
    amesos_ = rcp(new Amesos_Lapack((*reindexer_)(*lp_)));
  }

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
#if 0
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
        mllist.set("smoother: type "+(string)levelstr                    ,"MLS");
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
#endif
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
    case azprec_MLAPI:
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
    azlist.set("AZ_overlap",azvar->azoverlap);
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
    if (azvar->azprectyp == azprec_ML       ||
        azvar->azprectyp == azprec_MLfluid  ||
        azvar->azprectyp == azprec_MLfluid2 ||
        azvar->azprectyp == azprec_MLAPI       )
    {
      ParameterList& mllist = params.sublist("ML Parameters");
      ML_Epetra::SetDefaults("SA",mllist);
      switch (azvar->azprectyp)
      {
      case azprec_ML: // do nothing, this is standard
      break;
      case azprec_MLAPI: // set flag to use mlapi operator
        mllist.set<bool>("LINALG::AMG_Operator",true);
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
        {
          mllist.set("smoother: type "+(string)levelstr,"ILU");
          mllist.set("smoother: ifpack type "+(string)levelstr,"ILU");
          mllist.set("smoother: ifpack overlap "+(string)levelstr,azvar->azoverlap);
          mllist.set<double>("smoother: ifpack level-of-fill",(double)azvar->mlsmotimes[i]);
          ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
          ifpacklist.set<int>("fact: level-of-fill",azvar->mlsmotimes[i]);
          ifpacklist.set("schwarz: reordering type","rcm");
        }
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
        {
          mllist.set("coarse: type"          ,"IFPACK");
          mllist.set("coarse: ifpack type"   ,"ILU");
          mllist.set("coarse: ifpack overlap",0);
          ParameterList& ifpacklist = mllist.sublist("coarse: ifpack list");
          ifpacklist.set<int>("fact: level-of-fill",azvar->mlsmotimes[coarse]);
          ifpacklist.set("schwarz: reordering type","rcm");
        }
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



/*----------------------------------------------------------------------*
 | Multiply matrices A*B                                     mwgee 02/08|
 *----------------------------------------------------------------------*/
RCP<LINALG::SparseMatrix> LINALG::MLMultiply(const LINALG::SparseMatrix& A,
                                             const LINALG::SparseMatrix& B,
                                             bool complete)
{ 
  return MLMultiply(*A.EpetraMatrix(),*B.EpetraMatrix(),complete); 
}

/*----------------------------------------------------------------------*
 | Multiply matrices A*B                                     mwgee 02/08|
 *----------------------------------------------------------------------*/
//static void CopySortDeleteZeros(const Epetra_CrsMatrix& A, Epetra_CrsMatrix& As);
RCP<LINALG::SparseMatrix> LINALG::MLMultiply(const Epetra_CrsMatrix& A,
                                             const Epetra_CrsMatrix& B,
                                             bool complete)
{
  // make sure FillComplete was called on the matrices
  if (!A.Filled()) dserror("A has to be FillComplete");
  if (!B.Filled()) dserror("B has to be FillComplete");

  // For debugging, it might be helpful when all columns are
  // sorted and all zero values are wiped from the input:
  //RCP<Epetra_CrsMatrix> As = CreateMatrix(A.RowMap(),A.MaxNumEntries());
  //RCP<Epetra_CrsMatrix> Bs = CreateMatrix(B.RowMap(),B.MaxNumEntries());
  //CopySortDeleteZeros(A,*As);
  //CopySortDeleteZeros(B,*Bs);
  ML_Operator* ml_As = ML_Operator_Create(GetML_Comm());
  ML_Operator* ml_Bs = ML_Operator_Create(GetML_Comm());
  //ML_Operator_WrapEpetraMatrix(As.get(),ml_As);
  //ML_Operator_WrapEpetraMatrix(Bs.get(),ml_Bs);
  ML_Operator_WrapEpetraMatrix(const_cast<Epetra_CrsMatrix*>(&A),ml_As);
  ML_Operator_WrapEpetraMatrix(const_cast<Epetra_CrsMatrix*>(&B),ml_Bs);
  ML_Operator* ml_AtimesB = ML_Operator_Create(GetML_Comm());
  ML_2matmult(ml_As,ml_Bs,ml_AtimesB,ML_CSR_MATRIX);
  ML_Operator_Destroy(&ml_As);
  ML_Operator_Destroy(&ml_Bs);
  // For ml_AtimesB we have to reconstruct the column map in global indexing,
  // The following is going down to the salt-mines of ML ...
  int N_local = ml_AtimesB->invec_leng;
  ML_CommInfoOP* getrow_comm = ml_AtimesB->getrow->pre_comm;
  if (!getrow_comm) dserror("ML_Operator does not have CommInfo");
  ML_Comm* comm = ml_AtimesB->comm;
  if (N_local != B.DomainMap().NumMyElements()) 
    dserror("Mismatch in local row dimension between ML and Epetra");
  int N_rcvd  = 0;
  int N_send  = 0;
  int flag    = 0;
  for (int i=0; i<getrow_comm->N_neighbors; i++)
  {
    N_rcvd += (getrow_comm->neighbors)[i].N_rcv;
    N_send += (getrow_comm->neighbors)[i].N_send;
    if (  ((getrow_comm->neighbors)[i].N_rcv != 0) &&
       ((getrow_comm->neighbors)[i].rcv_list != NULL) )  flag = 1;
  }
  // For some unknown reason, ML likes to have stuff one larger than
  // neccessary...
  vector<double> dtemp(N_local+N_rcvd+1);
  vector<int> cmap(N_local+N_rcvd+1);
  for (int i=0; i<N_local; ++i)
  {
    cmap[i] = B.DomainMap().GID(i);
    dtemp[i] = (double)cmap[i];
  }
  ML_cheap_exchange_bdry(&dtemp[0],getrow_comm,N_local,N_send,comm);
  if (flag)
  {
    int count = N_local;
    const int neighbors = getrow_comm->N_neighbors;
    for (int i=0; i<neighbors; i++)
    {
      const int nrcv = getrow_comm->neighbors[i].N_rcv;
      for (int j=0; j<nrcv; j++)
        cmap[getrow_comm->neighbors[i].rcv_list[j]] = (int)dtemp[count++];
    }
  }
  else 
    for (int i=0; i<N_local+N_rcvd; ++i) cmap[i] = (int)dtemp[i];
  dtemp.clear();
  
  // we can now determine a matching column map for the result
  Epetra_Map gcmap(-1,N_local+N_rcvd,&cmap[0],0,A.Comm());

  // Allocate our result matrix and fill it
  // this is a very generous guess:
  int guessnpr = A.MaxNumEntries()*B.MaxNumEntries();
  RCP<Epetra_CrsMatrix> result 
    = rcp(new Epetra_CrsMatrix(Copy,A.RangeMap(),gcmap,guessnpr,false));

  int allocated=0;
  int rowlength;
  double* val=NULL;
  int* bindx=NULL;
  const int myrowlength = A.RowMap().NumMyElements();
  const Epetra_Map& rowmap = A.RowMap();
  vector<int> gcid(guessnpr);
  for (int i=0; i<myrowlength; ++i)
  {
    const int grid = rowmap.GID(i);
    // get local row
    ML_get_matrix_row(ml_AtimesB,1,&i,&allocated,&bindx,&val,&rowlength,0);
    if (!rowlength) continue;
    if ((int)gcid.size() < rowlength) gcid.resize(rowlength);
    for (int j=0; j<rowlength; ++j)
    {
      gcid[j] = gcmap.GID(bindx[j]);
#ifdef DEBUG
      if (gcid[j]<0) dserror("This is really bad... cannot find gcid");
#endif
    }
#ifdef DEBUG
    int err = result->InsertGlobalValues(grid,rowlength,val,&gcid[0]);
    if (err!=0 && err!=1) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
#else
    result->InsertGlobalValues(grid,rowlength,val,&gcid[0]);
#endif
  }
  if (bindx) ML_free(bindx);
  if (val) ML_free(val);
  ML_Operator_Destroy(&ml_AtimesB);
  if (complete) 
  {
    int err = result->FillComplete(B.DomainMap(),A.RangeMap());
    if (err) dserror("Epetra_CrsMatrix::FillComplete returned err=%d",err);
  }
  return rcp(new SparseMatrix(result));
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/*
static void CopySortDeleteZeros(const Epetra_CrsMatrix& A, Epetra_CrsMatrix& As)
{
  vector<int>    scindices(A.MaxNumEntries());
  vector<double> scvalues(A.MaxNumEntries());
  for (int i=0; i<A.NumMyRows(); ++i)
  {
    int grid = A.RowMap().GID(i);
    int numentries;
    double* values;
    int* indices;
    A.ExtractMyRowView(i,numentries,values,indices);
    int snumentries=0;
    for (int j=0; j<numentries; ++j)
    {
      if (values[j]==0.0) continue;
      scindices[snumentries] = A.ColMap().GID(indices[j]);
      scvalues[snumentries] = values[j];
      snumentries++;
    }
    ML_az_sort(&scindices[0],snumentries,NULL,&scvalues[0]);
    int err = As.InsertGlobalValues(grid,snumentries,&scvalues[0],&scindices[0]);
    if (err) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
  }
  if (A.Filled()) As.FillComplete(A.DomainMap(),A.RangeMap(),true);
  return;
}
*/


#endif  // #ifdef CCADISCRET
