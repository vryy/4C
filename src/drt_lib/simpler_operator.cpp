/*!----------------------------------------------------------------------
\file simpler_operator.cpp

\class LINALG::SIMPLER_Operator

\brief An approximate block factorization preconditioner based on the
       SIMPLE family of methods

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#include "simpler_operator.H"

#define SIMPLEC_DIAGONAL      0    // 1: row sums     0: just diagonal
#define CHEAPSIMPLE_ALGORITHM 1    // 1: AMG          0: true solve
#define SIMPLER_ALGORITHM     0    // 1: triple solve 0: double solve
#define SIMPLER_ALPHA         1.0  // simple pressure damping parameter
#define SIMPLER_TIMING        0    // printout timing of setup
/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 02/08|
 *----------------------------------------------------------------------*/
LINALG::SIMPLER_Operator::SIMPLER_Operator(RCP<Epetra_CrsMatrix> A,
                                           const ParameterList& velocitylist,
                                           const ParameterList& pressurelist,
                                           FILE* outfile)
  : outfile_(outfile),
    vlist_(velocitylist),
    plist_(pressurelist),
    alpha_(SIMPLER_ALPHA)
{
  // remove the SIMPLER sublist from the vlist_,
  // otherwise it will try to recursively create a SIMPLE
  // preconditioner when we do the subblock solvers
  if (vlist_.isSublist("SIMPLER")) vlist_.remove("SIMPLER");

  Setup(A,velocitylist,pressurelist);

  return;
}


/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::SIMPLER_Operator::Setup(RCP<Epetra_CrsMatrix> A,
                                     const ParameterList& origvlist,
                                     const ParameterList& origplist)
{
  const int myrank = A->Comm().MyPID();
  Epetra_Time time(A->Comm());
  Epetra_Time totaltime(A->Comm());
  
  // see whether velocity and pressure solver where configured as ML
  bool visml = vlist_.isSublist("ML Parameters");
  bool pisml = plist_.isSublist("ML Parameters");
  if (!visml) dserror("SIMPLER only works with ML-AMG for velocity");

  // get # dofs per node from vlist_ and split row map
  time.ResetStartTime();
  const int ndofpernode = vlist_.sublist("ML Parameters").get<int>("PDE equations",0);
  if (ndofpernode != 4 && ndofpernode !=3) dserror("You should have either 3 or 4 dofs per node at this point");
  const Epetra_Map& fullmap = A->RowMap();
  const int length = fullmap.NumMyElements();
  const int nv     = ndofpernode-1;
  const int nlnode = length / ndofpernode;
  vector<int> vgid(nlnode*nv);
  vector<int> pgid(nlnode);
  int vcount=0;
  for (int i=0; i<nlnode; ++i)
  {
    for (int j=0; j<ndofpernode-1; ++j)
      vgid[vcount++] = fullmap.GID(i*ndofpernode+j);
    pgid[i] = fullmap.GID(i*ndofpernode+ndofpernode-1);
  }
  vector<RCP<const Epetra_Map> > maps(2);
  maps[0] = rcp(new Epetra_Map(-1,nlnode*nv,&vgid[0],0,fullmap.Comm()));
  maps[1] = rcp(new Epetra_Map(-1,nlnode,&pgid[0],0,fullmap.Comm()));
  vgid.clear(); pgid.clear();
  mmex_.Setup(fullmap,maps);
  if (!myrank && SIMPLER_TIMING) printf("--- Time to split map       %10.3E\n",time.ElapsedTime());
  time.ResetStartTime();

  // wrap matrix in SparseMatrix and split it into 2x2 BlockMatrix
  {
    SparseMatrix fullmatrix(A);
    A_ = fullmatrix.Split<LINALG::DefaultBlockMatrixStrategy>(mmex_,mmex_);
    if (!myrank && SIMPLER_TIMING) printf("--- Time to split matrix    %10.3E\n",time.ElapsedTime());
    time.ResetStartTime();
    A_->Complete();
    if (!myrank && SIMPLER_TIMING) printf("--- Time to complete matrix %10.3E\n",time.ElapsedTime());
    time.ResetStartTime();
  }

  // split nullspace into velocity and pressure subproblem
  if (visml)
  {
    vlist_.sublist("ML Parameters").set("PDE equations",nv);
    vlist_.sublist("ML Parameters").set("null space: dimension",nv);
    const int vlength = (*A_)(0,0).RowMap().NumMyElements();
    RCP<vector<double> > vnewns = rcp(new vector<double>(nv*vlength,0.0));
    for (int i=0; i<nlnode; ++i)
    {
      (*vnewns)[i*nv] = 1.0;
      (*vnewns)[vlength+i*nv+1] = 1.0;
      (*vnewns)[2*vlength+i*nv+2] = 1.0;
    }
    vlist_.sublist("ML Parameters").set("null space: vectors",&((*vnewns)[0]));
    vlist_.sublist("ML Parameters").remove("nullspace",false);
    vlist_.sublist("Michael's secret vault").set<RCP<vector<double> > >("velocity nullspace",vnewns);
  }
  if (!myrank && SIMPLER_TIMING) printf("--- Time to do v nullspace  %10.3E\n",time.ElapsedTime());
  time.ResetStartTime();

  if (pisml)
  {
    plist_.sublist("ML Parameters").set("PDE equations",1);
    plist_.sublist("ML Parameters").set("null space: dimension",1);
    const int plength = (*A_)(1,1).RowMap().NumMyElements();
    RCP<vector<double> > pnewns = rcp(new vector<double>(plength,1.0));
    plist_.sublist("ML Parameters").set("null space: vectors",&((*pnewns)[0]));
    plist_.sublist("ML Parameters").remove("nullspace",false);
    plist_.sublist("Michael's secret vault").set<RCP<vector<double> > >("pressure nullspace",pnewns);
  }
  if (!myrank && SIMPLER_TIMING) printf("--- Time to do p nullspace  %10.3E\n",time.ElapsedTime());
  time.ResetStartTime();

  // Modify lists to reuse subblock preconditioner at least maxiter times
  {
    int maxiter = vlist_.sublist("Aztec Parameters").get("AZ_max_iter",1);
    vlist_.sublist("Aztec Parameters").set("reuse",maxiter+1);
    plist_.sublist("Aztec Parameters").set("reuse",maxiter+1);
  }

#if SIMPLEC_DIAGONAL
  // Allocate and compute abs(rowsum(A(0,0))^{-1}
  {
    Epetra_Vector diag(*mmex_.Map(0),true);
    RCP<Epetra_CrsMatrix> A00 = (*A_)(0,0).EpetraMatrix();
    const int nlrows = A00->RowMap().NumMyElements();
    for (int i=0; i<nlrows; ++i)
    {
      int numentries;
      double* vals;
      int err = A00->ExtractMyRowView(i,numentries,vals);
      if (err) dserror("Epetra_CrsMatrix::ExtractMyRowView returned %d",err);
      for (int j=0; j<numentries; ++j) diag[i] += vals[j];
      diag[i] = abs(diag[i]);
    }
    int err = diag.Reciprocal(diag);
    if (err) dserror("Epetra_MultiVector::Reciprocal returned %d",err);
    diagAinv_ = rcp(new SparseMatrix(diag));
    diagAinv_->Complete(*mmex_.Map(0),*mmex_.Map(0));
  }
#else
  // Allocate and compute diag(A(0,0)^{-1}
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

  // Allocate and compute approximate Schur complement operator S
  // S = A(1,1) - A(1,0) * diagAinv * A(0,1)
  {
    Epetra_Time ltime(A_->Comm());
    S_ = LINALG::Multiply(*diagAinv_,false,(*A_)(0,1),false,true);
    //cout << *S_; exit(0);
    if (!myrank && SIMPLER_TIMING) printf("*** S = diagAinv * A(0,1) %10.3E\n",ltime.ElapsedTime());
    ltime.ResetStartTime();
    S_ = LINALG::MLMultiply((*A_)(1,0),*S_,false);
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
  // Allocate preconditioner for pressure and velocity
  {
    if (!visml || !pisml) dserror("Have to use ML for this variant of SIMPLE");
    Pv_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*((*A_)(0,0).EpetraMatrix()),
                                                      vlist_.sublist("ML Parameters"),
                                                      true));
    if (!myrank && SIMPLER_TIMING) printf("--- Time to do AMG(v)       %10.3E\n",time.ElapsedTime());
    time.ResetStartTime();
    Pp_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*(S_->EpetraMatrix()),
                                                      plist_.sublist("ML Parameters"),
                                                      true));
    if (!myrank && SIMPLER_TIMING) printf("--- Time to do AMG(p)       %10.3E\n",time.ElapsedTime());
    time.ResetStartTime();
  }
#else
  // Allocate solver for pressure and velocity
  {
    RCP<ParameterList> vrcplist = rcp(&vlist_,false);
    vsolver_ = rcp(new LINALG::Solver(vrcplist,A_->Comm(),outfile_));
    RCP<ParameterList> prcplist = rcp(&plist_,false);
    psolver_ = rcp(new LINALG::Solver(prcplist,A_->Comm(),outfile_));
  }
#endif

  // Allocate velocity and pressure solution and rhs vectors
  vx_ = LINALG::CreateVector(*mmex_.Map(0),false);
  vb_ = LINALG::CreateVector(*mmex_.Map(0),false);
  px_ = LINALG::CreateVector(*mmex_.Map(1),false);
  pb_ = LINALG::CreateVector(*mmex_.Map(1),false);

  // Allocate working vectors for velocity and pressure
  vwork1_  = LINALG::CreateVector(*mmex_.Map(0),false);
  vwork2_  = LINALG::CreateVector(*mmex_.Map(0),false);
  pwork1_ = LINALG::CreateVector(*mmex_.Map(1),false);
  pwork2_ = LINALG::CreateVector(*mmex_.Map(1),false);

  if (!myrank && SIMPLER_TIMING) printf("--- Time to do allocate mem %10.3E\n",time.ElapsedTime());
  if (!myrank && SIMPLER_TIMING) printf("=== Total simpler setup === %10.3E\n",totaltime.ElapsedTime());

  return;
}


/*----------------------------------------------------------------------*
 |  apply const operator (public)                            mwgee 02/08|
 *----------------------------------------------------------------------*/
int LINALG::SIMPLER_Operator::ApplyInverse(const Epetra_MultiVector& X,
                                           Epetra_MultiVector& Y) const
{
  // note: Aztec might pass X and Y as physically identical objects,
  // so we better deep copy here

  // extract initial guess and rhs for velocity and pressure
  mmex_.ExtractVector(X,0,*vb_);
  mmex_.ExtractVector(X,1,*pb_);

#if CHEAPSIMPLE_ALGORITHM // SIMPLE and SIMPLEC but without solve, just AMG
  CheapSimple(vx_,px_,vb_,pb_);
#else

#if SIMPLER_ALGORITHM
  Simpler(vx_,px_,vb_,pb_);
#else
  Simple(vx_,px_,vb_,pb_);
#endif

#endif

  // insert solution for velocity and pressure
  mmex_.InsertVector(*vx_,0,Y);
  mmex_.InsertVector(*px_,1,Y);

  return 0;
}



/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::SIMPLER_Operator::Simple(RCP<Epetra_Vector> vx, RCP<Epetra_Vector> px,
                                      RCP<Epetra_Vector> vb, RCP<Epetra_Vector> pb) const
{

  //------------------------------------------------------------ L-solve
  // Solve A(0,0) \hat{v} = vb , result is on vwork1_
  vwork1_->PutScalar(0.0);
  vsolver_->Solve((*A_)(0,0).EpetraMatrix(),vwork1_,vb,true,false);

  // Build rhs for second solve pwork1_ = pb - A(1,0)*\hat{v} , result is on pwork1_
  (*A_)(1,0).Multiply(false,*vwork1_,*pwork1_);
  pwork1_->Update(1.0,*pb,-1.0);

  // Solve S \hat{p} = pb - A(1,0)*\hat{v} , result is on px
  px->PutScalar(0.0);
  psolver_->Solve(S_->EpetraMatrix(),px,pwork1_,true,false);

  //------------------------------------------------------------ U-solve
  // p = 1/alpha * \hat{p} , result is on px
  if (alpha_ != 1.0) px->Scale(1./alpha_);

  // v = \hat{v} - diag(A(0,0))^{-1} A(0,1) p , result is on vx
  (*A_)(0,1).Multiply(false,*px,*vwork2_);
  diagAinv_->Multiply(false,*vwork2_,*vx);
  vx->Update(1.0,*vwork1_,-1.0);

  return;
}



/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::SIMPLER_Operator::Simpler(RCP<Epetra_Vector> vx, RCP<Epetra_Vector> px,
                                       RCP<Epetra_Vector> vb, RCP<Epetra_Vector> pb) const
{

  //------------------------------------------------------------ L-solve
  // Solve I \hat{v} = vb, result is on vb
  // nothing to do

  // Build rhs for second solve pwork1_ = pb - A(1,0) diag(A(0,0))^{-1} \hat{v}
  diagAinv_->Multiply(false,*vb,*vwork1_);
  (*A_)(1,0).Multiply(false,*vwork1_,*pwork1_);
  pwork1_->Update(1.0,*pb,-1.0);

  // Solve I \hat{p} = pb - A(1,0) diag(A(0,0))^{-1} \hat{v}, result is on pwork1_
  // nothing to do

  //------------------------------------------------------------ U-solve
  // Solve S p = \hat{p}, result is on px
  px->PutScalar(0.0);
  psolver_->Solve(S_->EpetraMatrix(),px,pwork1_,true,false);

  // Build rhs for second solve \hat{v} - A(0,1) p, result is on vb
  (*A_)(0,1).Multiply(false,*px,*vwork1_);
  vb->Update(-1.0,*vwork1_,1.0);

  // Solve A(0,0) v = \hat{v} - A(0,1) p, result is on vx
  vx->PutScalar(0.0);
  vsolver_->Solve((*A_)(0,0).EpetraMatrix(),vx,vb,true,false);

  //------------------------------------------------ Implicit projection
  // px = \alpha I px
  if (alpha_ != 1.0) px->Scale(alpha_);

  // vx = ( I + diag(A(0,0))^{-1} A(0,1) S^-1 A(1,0) ) vx
  (*A_)(1,0).Multiply(false,*vx,*pwork1_);
  pwork2_->PutScalar(0.0);
  psolver_->Solve(S_->EpetraMatrix(),pwork2_,pwork1_,false,false);
  (*A_)(0,1).Multiply(false,*pwork2_,*vwork1_);
  diagAinv_->Multiply(false,*vwork1_,*vwork2_);
  vx->Update(-1.0,*vwork2_,1.0);

  return;
}


/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 02/08|
 *----------------------------------------------------------------------*/
void LINALG::SIMPLER_Operator::CheapSimple(RCP<Epetra_Vector> vx, RCP<Epetra_Vector> px,
                                           RCP<Epetra_Vector> vb, RCP<Epetra_Vector> pb) const
{
  //------------------------------------------------------------ L-solve
  // Solve A(0,0) \hat{v} = vb , result is on vwork1_
  Pv_->ApplyInverse(*vb,*vwork1_);

  // Build rhs for second solve pwork1_ = pb - A(1,0)*\hat{v} , result is on pwork1_
  (*A_)(1,0).Multiply(false,*vwork1_,*pwork1_);
  pwork1_->Update(1.0,*pb,-1.0);

  // Solve S \hat{p} = pb - A(1,0)*\hat{v} , result is on px
  Pp_->ApplyInverse(*pwork1_,*px);

  //------------------------------------------------------------ U-solve
  // p = 1/alpha * \hat{p} , result is on px
  if (alpha_ != 1.0) px->Scale(1./alpha_);

  // v = \hat{v} - diag(A(0,0))^{-1} A(0,1) p , result is on vx
  (*A_)(0,1).Multiply(false,*px,*vwork2_);
  diagAinv_->Multiply(false,*vwork2_,*vx);
  vx->Update(1.0,*vwork1_,-1.0);

  return;
}









#endif  // #ifdef CCADISCRET
