#include "vm3_solver.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 03/06|
 *----------------------------------------------------------------------*/
VM3_Solver::VM3_Solver(
                                    RefCountPtr<Epetra_CrsMatrix> Aplus,
                                    RefCountPtr<Epetra_CrsMatrix> A,
                                    ParameterList& mlparams,
                                    bool compute) :
iscomputed_(false),
mlparams_(mlparams),
Aplus_(Aplus),
A_(A)
{
  label_  = "VM3_Solver";


  if (compute) Compute();
    
  return;
}

/*----------------------------------------------------------------------*
 |  FAS multigrid solver for VM3                             m.gee 03/06|
 *----------------------------------------------------------------------*/
int VM3_Solver::Solve(
                     const Epetra_MultiVector& B, Epetra_MultiVector& X)
{
    // do setup if not already done
  if (!iscomputed_)
  {
    VM3_Solver& tmp = 
                       const_cast<VM3_Solver&>(*this);
    tmp.Compute();
  }


  // create a Space
  const Epetra_BlockMap& bmap = B.Map();
  Space space;
  space.Reshape(bmap.NumGlobalElements(),bmap.NumMyElements(),bmap.MyGlobalElements());
  
  // create input/output mlapi multivectors
  MultiVector b_f(space,1,false);
  MultiVector x_f(space,1,false);
  const int nele = B.Map().NumMyElements();
  for (int i=0; i<nele; ++i)
  {
    x_f(i) = X[0][i];
    b_f(i) = B[0][i];
  }

  // call AMG
  MultiLevelVCycle(b_f,x_f);

  // copy solution back
  for (int i=0; i<nele; ++i)
    X[0][i] = x_f(i);
  
  return 0;

}

/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (public)           m.gee 03/06|
 *----------------------------------------------------------------------*/
int VM3_Solver::ApplyInverse(
                     const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // apply the preconditioner to X and return result in Y
  
  // do setup if not already done
  if (!iscomputed_)
  {
    VM3_Solver& tmp = 
                       const_cast<VM3_Solver&>(*this);
    tmp.Compute();
  }


  // create a Space
  const Epetra_BlockMap& bmap = X.Map();
  Space space;
  space.Reshape(bmap.NumGlobalElements(),bmap.NumMyElements(),bmap.MyGlobalElements());
  
  // create input/output mlapi multivectors
  MultiVector b_f(space,1,false);
  MultiVector x_f(space,1,false);
  const int nele = X.Map().NumMyElements();
  for (int i=0; i<nele; ++i)
  {
    x_f(i) = Y[0][i];
    b_f(i) = X[0][i];
  }

  // call AMG
  MultiLevelVCycle(b_f,x_f);

  // copy solution back
  for (int i=0; i<nele; ++i)
    Y[0][i] = x_f(i);
  
  return 0;
}

/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (private)          m.gee 03/06|
 *----------------------------------------------------------------------*/
int VM3_Solver::MultiLevelVCycle(MultiVector& b_f,
                                 MultiVector& x_f)
const
{
  int level, levelm1;

  // smoothing on finest level
  {
      // multivector definitions
      MultiVector x_p0(P(0).GetRangeSpace(),1,false);

      // step 1: scale separation
      // scale part on level 1 prolongated to level 0
      x_p0 = P(0) * R(0) * x_f;

      // scale part on current level 0
      x_f = x_f - x_p0;

      // step 2: RHS computation
      // compute additional RHS-term for scale part on level 1
      b_f = b_f - A(0) * x_p0;

      // step 3: pre-smoothing
      S(0).Apply(b_f,x_f);

      // step 4: composition of complete solution
      x_f = x_p0 + x_f;
  }

  // smoothing on medium levels
  for (level=1; level<maxlevels_-2; ++level)
  {
      levelm1=level-1;

      // multivector definitions
      MultiVector b_r(P(level).GetRangeSpace(),1,false);
      MultiVector x_r(P(level).GetRangeSpace(),1,false);
      MultiVector x_p(P(level).GetRangeSpace(),1,false);
      MultiVector b_c(P(level).GetRangeSpace(),1,false);
      MultiVector x_c(P(level).GetRangeSpace(),1,false);
      MultiVector x_rm(P(levelm1).GetRangeSpace(),1,false);
      MultiVector x_m(P(levelm1).GetRangeSpace(),1,false);
      MultiVector x_p0(P(0).GetRangeSpace(),1,false);
      MultiVector x_c0(P(0).GetRangeSpace(),1,false);
      MultiVector x_m0(P(0).GetRangeSpace(),1,false);

      // pre-step: solution and RHS restricted to current level
      Restrict(x_f,x_r,level);

      // step 1: scale separation
      // scale part on level+1 prolongated to current level
      x_p = P(level) * R(level) * x_r;
      // and prolongated to finest level
      Prolong(x_p,x_p0,level);

      // scale part on current level
      x_c = x_r - x_p;

      // (combined) scale parts on finer levels restricted to level-1
      Restrict(x_f,x_rm,levelm1);
      x_m = x_rm - P(levelm1) * R(levelm1) * x_rm;
      // and prolongated to finest level if necessary
      if (level > 1)
      {
         Prolong(x_m,x_m0,levelm1);
      }

      // step 2: RHS computation
      // compute additional RHS-terms for scale parts on level+1 and finer levels
      Restrict(b_f,b_r,level);
      b_c = b_r - A(level) * x_p - R(levelm1) * A(levelm1) * x_m;

      // step 3: pre-smoothing
      S(level).Apply(b_c,x_c);

      // step 4: composition of complete solution
      // before: result for current scale part prolongated to finest level
      Prolong(x_c,x_c0,level);
      x_f = x_p0 + x_c0 + x_m0;
  }

  // solution on coarsest level
  {
      level = maxlevels_-1;
      levelm1=level-1;

      // multivector definitions
      MultiVector b_r(P(level).GetRangeSpace(),1,false);
      MultiVector x_r(P(level).GetRangeSpace(),1,false);
      MultiVector b_c(P(level).GetRangeSpace(),1,false);
      MultiVector x_rm(P(levelm1).GetRangeSpace(),1,false);
      MultiVector x_m(P(levelm1).GetRangeSpace(),1,false);
      MultiVector x_c0(P(0).GetRangeSpace(),1,false);
      MultiVector x_m0(P(0).GetRangeSpace(),1,false);

      // pre-step: solution and RHS restricted to current level
      Restrict(x_f,x_r,level);

      // step 1: scale separation
      // (combined) scale parts on finer levels restricted to level-1
      Restrict(x_f,x_rm,levelm1);
      x_m = x_rm - P(levelm1) * R(levelm1) * x_rm;
      // and prolongated to finest level if necessary
      if (level > 1)
      {
         Prolong(x_m,x_m0,levelm1);
      }

      // step 2: RHS computation
      // compute additional RHS-terms for scale parts on finer levels
      Restrict(b_f,b_r,level);
      b_c = b_r - R(levelm1) * A(levelm1) * x_m;

      // step 3: solution
      x_r = S(level) * b_c;

      // step 4: composition of complete solution
      // before: result for current scale part prolongated to finest level
      Prolong(x_r,x_c0,level);
      x_f = x_c0 + x_m0;
  }

  // solution on medium levels
  for (level=maxlevels_-2; level>0; --level)
  {
      levelm1=level-1;

      // multivector definitions
      MultiVector b_r(P(level).GetRangeSpace(),1,false);
      MultiVector x_r(P(level).GetRangeSpace(),1,false);
      MultiVector x_p(P(level).GetRangeSpace(),1,false);
      MultiVector b_c(P(level).GetRangeSpace(),1,false);
      MultiVector x_c(P(level).GetRangeSpace(),1,false);
      MultiVector x_rm(P(levelm1).GetRangeSpace(),1,false);
      MultiVector x_m(P(levelm1).GetRangeSpace(),1,false);
      MultiVector x_p0(P(0).GetRangeSpace(),1,false);
      MultiVector x_c0(P(0).GetRangeSpace(),1,false);
      MultiVector x_m0(P(0).GetRangeSpace(),1,false);

      // pre-step: solution and RHS restricted to current level
      Restrict(x_f,x_r,level);

      // step 1: scale separation
      // scale part on level+1 prolongated to current level
      x_p = P(level) * R(level) * x_r;
      // and prolongated to finest level
      Prolong(x_p,x_p0,level);

      // scale part on current level
      x_c = x_r - x_p;

      // (combined) scale parts on finer levels restricted to level-1
      Restrict(x_f,x_rm,levelm1);
      x_m = x_rm - P(levelm1) * R(levelm1) * x_rm;
      // and prolongated to finest level if necessary
      if (level > 2)
      {
         Prolong(x_m,x_m0,levelm1);
      }

      // step 2: RHS computation
      // compute additional RHS-terms for scales on level+1 and finer levels
      Restrict(b_f,b_r,level);
      b_c = b_r - A(level) * x_p - R(levelm1) * A(levelm1) * x_m;

      // step 3: pre-smoothing
      S(level).Apply(b_c,x_c);

      // step 4: composition of complete solution
      // before: result for current scale part prolongated to finest level
      Prolong(x_c,x_c0,level);
      x_f = x_p0 + x_c0 + x_m0;
  }

  // smoothing on finest level
  {
      // multivector definitions
      MultiVector x_p0(P(0).GetRangeSpace(),1,false);

      // step 1: scale separation
      // scale part on level 1 prolongated to level 0
      x_p0 = P(0) * R(0) * x_f;

      // scale part on current level 0
      x_f = x_f - x_p0;

      // step 2: RHS computation
      // compute additional RHS-term for scale part on level 1
      b_f = b_f - A(0) * x_p0;

      // step 3: pre-smoothing
      S(0).Apply(b_f,x_f);

      // step 4: composition of complete solution
      x_f = x_p0 + x_f;
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  restriction operation (private)                             vg 08/07|
 *----------------------------------------------------------------------*/
int VM3_Solver::Restrict(const MultiVector& x_f,
                         MultiVector& x_r,
                         int   level)
const
{
  MultiVector x_rf(P(0).GetRangeSpace(),1,false);
  MultiVector x_rc(P(level).GetRangeSpace(),1,false);

  x_rf = x_f;

  for (int i=0; i<level; ++i)
  {
      MultiVector x_rf(P(level).GetRangeSpace(),1,false);
      MultiVector x_rc(P(level).GetDomainSpace(),1,false);

      x_rc = R(i) * x_rf;
  }

  x_r = x_rc;

  return 0;
}

/*----------------------------------------------------------------------*
 |  prolongation operation (private)                            vg 08/07|
 *----------------------------------------------------------------------*/
int VM3_Solver::Prolong(const MultiVector& x_l,
                        MultiVector& x_0,
                        int   level)
const
{
  MultiVector x_pc(P(level).GetRangeSpace(),1,false);
  MultiVector x_pf(P(0).GetRangeSpace(),1,false);

  x_pc = x_l;

  for (int i=0; i<level; ++i)
  {
      MultiVector x_pf(P(level-1-i).GetRangeSpace(),1,false);
      MultiVector x_pc(P(level-1-i).GetDomainSpace(),1,false);

      x_pf = P(level-1-i) * x_pc;
  }

  x_0 = x_pf;

  return 0;
}

/*----------------------------------------------------------------------*
 |  compute the preconditioner (public)                      m.gee 03/06|
 *----------------------------------------------------------------------*/
bool VM3_Solver::Compute()
{
  // setup phase of multigrid
  iscomputed_ = false;
  
  // this is important to have!!!
  MLAPI::Init();
  
  // get parameters
  int     maxlevels     = mlparams_.get("max levels",10);
  int     maxcoarsesize = mlparams_.get("coarse: max size",10);
  double* nullspace     = mlparams_.get("null space: vectors",(double*)NULL);
  int     nsdim         = mlparams_.get("null space: dimension",1);
  int     numpde        = mlparams_.get("PDE equations",1);
  double  damping       = mlparams_.get("aggregation: damping factor",1.33);
  string  eigenanalysis = mlparams_.get("eigen-analysis: type", "Anorm");
  string  ptype         = mlparams_.get("prolongator: type","mod_full");

  // Enhanced parameter extraction in case of list generated via ML-Prec.
  // string  fsmoothertype = mlparams_.get("smoother: type (level 0)","symmetric Gauss-Seidel");
  // double  fsmootherdamp = mlparams_.get("smoother: damping factor (level 0)",1.33);
  // int     fsmoothsweeps = mlparams_.get("smoother: sweeps (level 0)",2);
  // string  smoothertype  = mlparams_.get("smoother: type (level 1)","symmetric Gauss-Seidel");
  // double  smootherdamp  = mlparams_.get("smoother: damping factor (level 1)",1.33);
  // int     smoothsweeps  = mlparams_.get("smoother: sweeps (level 1)",2);
  string  coarsetype    = mlparams_.get("coarse: type","Amesos-KLU");
  // double  coarsedamp    = mlparams_.get("coarse: damping factor",1.33);
  // int     coarsesweeps  = mlparams_.get("coarse: sweeps",2);

  // Currently, only one smoother type can be extracted anyway (in MLAPI::InverseOperator::Reshape), 
  // so it doesn't make sense to include different relaxations here. Similarly, only
  // one smoother type is included in the parameter list mlparams. See linalg_solver.cpp
  string  smoothertype  = mlparams_.get("smoother: type","symmetric Gauss-Seidel");
  string  fsmoothertype = smoothertype;

  Space space(A_->RowMatrixRowMap());
  Operator mlapiA(space,space,A_.get(),false);
  Operator mlapiAplus(space,space,Aplus_.get(),false);


  mlapiRmod_.resize(maxlevels);                       
  mlapiPmod_.resize(maxlevels);                       
  //mlapiRP_.resize(maxlevels);                       
  //mlapiPR_.resize(maxlevels);
  mlapiRA_.resize(maxlevels);
  mlapiA_.resize(maxlevels);
  mlapiS_.resize(maxlevels);
  mlapiAplus_.resize(1);

  // build nullspace;
  MultiVector NS;
  MultiVector NextNS;
  
  NS.Reshape(mlapiA.GetRangeSpace(),nsdim);
  if (nullspace)
  {
    for (int i=0; i<nsdim; ++i)
      for (int j=0; j<NS.GetMyLength(); ++j)
        NS(j,i) = nullspace[i*NS.GetMyLength()+j];
  }
  else
  {
    if (numpde==1) NS = 1.0;
    else
    {
      NS = 0.0;
      for (int i=0; i<NS.GetMyLength(); ++i)
        for (int j=0; j<numpde; ++j)
          if ( i % numpde == j)
            NS(i,j) = 1.0;
    }
  }

  double lambdamax;
  Operator Ptent;
  Operator P;
  Operator Rtent;
  Operator R;
  Operator IminusA;
  Operator C;

  Operator Pmod;
  Operator Rmod;
  InverseOperator S;

  mlapiAplus_[0] = mlapiAplus;
  mlapiA_[0] = mlapiA;

  // build the operators for level 0 first
  int level = 0;

  // build smoother
  if (Comm().MyPID()==0)
  {
    ML_print_line("-", 78);
    cout << "VM3 solver : creating smoother level " << level << endl; 
    fflush(stdout);
  }
  S.Reshape(mlapiAplus,fsmoothertype,mlparams_);
  
  if (level) mlparams_.set("PDE equations", NS.GetNumVectors());
  
 
  if (Comm().MyPID()==0)
  {
    ML_print_line("-", 80);
    cout << "VM3 solver: creating level " << level+1 << endl;
    ML_print_line("-", 80);
    fflush(stdout);
  }
  mlparams_.set("workspace: current level",level);
  GetPtent(mlapiA,mlparams_,NS,Ptent,NextNS);
  NS = NextNS;
  
  if (damping)
  {
    if (eigenanalysis == "Anorm")
      lambdamax = MaxEigAnorm(mlapiA,true);
    else if (eigenanalysis == "cg")
      lambdamax = MaxEigCG(mlapiA,true);
    else if (eigenanalysis == "power-method")
      lambdamax = MaxEigPowerMethod(mlapiA,true);
    else ML_THROW("incorrect parameter (" + eigenanalysis + ")", -1);
  
    IminusA = GetJacobiIterationOperator(mlapiA,damping/lambdamax);
    P = IminusA * Ptent;
  }
  else
  {
    P = Ptent;
    lambdamax = -1.0;
  }
  
  R = GetTranspose(P);
  if (damping)
    Rtent = GetTranspose(Ptent);
  else
    Rtent = R;
    
  // variational coarse grid
  C = GetRAP(R,mlapiA,P);

  // build the matrix-matrix products R*A, R*P and P*R
  ML_Operator* Rmat = R.GetML_Operator();
  ML_Operator* Amat = mlapiA.GetML_Operator();
  ML_Operator* Pmat = P.GetML_Operator();
  ML_Operator* RAmat;
//  ML_Operator* RPmat;
//  ML_Operator* PRmat;
  ML_matmat_mult(Rmat, Amat, &RAmat);
//  ML_matmat_mult(Rmat, Pmat, &RPmat);
//  ML_matmat_mult(Pmat, Rmat, &PRmat);
  Operator RA(mlapiA.GetDomainSpace(),R.GetRangeSpace(), RAmat,false);
//  Operator RP(P.GetDomainSpace(),R.GetRangeSpace(), RPmat,false);
//  Operator PR(R.GetDomainSpace(),P.GetRangeSpace(), PRmat,false);
  
  // write the temporary values into the correct slot
  mlapiRA_[level]       = RA;
//  mlapiRP_[level]       = RP;
//  mlapiPR_[level]       = PR;
  mlapiRmod_[level]     = R;  
  mlapiPmod_[level]     = P;
  mlapiA_[level+1]      = C;
  mlapiS_[level]        = S;

  // loop for level 1 -- maxleves-1
  for (level=1; level<maxlevels-1; ++level)
  {
    // this level's operator
    mlapiA = mlapiA_[level];

    // build smoother
    if (Comm().MyPID()==0)
    {
      ML_print_line("-", 78);
      cout << "VM3 solver : creating smoother level " << level << endl; 
      fflush(stdout);
    }
    S.Reshape(mlapiA,smoothertype,mlparams_);
    
    if (level) mlparams_.set("PDE equations", NS.GetNumVectors());
    
  
    if (Comm().MyPID()==0)
    {
      ML_print_line("-", 80);
      cout << "VM3 solver : creating level " << level+1 << endl;
      ML_print_line("-", 80);
      fflush(stdout);
    }

    mlparams_.set("workspace: current level",level);
    GetPtent(mlapiA,mlparams_,NS,Ptent,NextNS);
    NS = NextNS;
    
    if (damping)
    {
      if (eigenanalysis == "Anorm")
        lambdamax = MaxEigAnorm(mlapiA,true);
      else if (eigenanalysis == "cg")
        lambdamax = MaxEigCG(mlapiA,true);
      else if (eigenanalysis == "power-method")
        lambdamax = MaxEigPowerMethod(mlapiA,true);
      else ML_THROW("incorrect parameter (" + eigenanalysis + ")", -1);
    
      IminusA = GetJacobiIterationOperator(mlapiA,damping/lambdamax);
      P = IminusA * Ptent;
    }
    else
    {
      P = Ptent;
      lambdamax = -1.0;
    }
    
    R = GetTranspose(P);
    if (damping)
      Rtent = GetTranspose(Ptent);
    else
      Rtent = R;
      
    // variational coarse grid
    C = GetRAP(R,mlapiA,P);

    // build the matrix-matrix products R*A, R*P and P*R
    Rmat = R.GetML_Operator();
    Amat = mlapiA.GetML_Operator();
    Pmat = P.GetML_Operator();
    ML_matmat_mult(Rmat, Amat, &RAmat);
//    ML_matmat_mult(Rmat, Pmat, &RPmat);
//    ML_matmat_mult(Pmat, Rmat, &PRmat);
    Operator RA(mlapiA.GetDomainSpace(),R.GetRangeSpace(), RAmat,false);
//    Operator RP(P.GetDomainSpace(),R.GetRangeSpace(), RPmat,false);
//    Operator PR(R.GetDomainSpace(),P.GetRangeSpace(), PRmat,false);

    // write the temporary values into the correct slot
    mlapiRA_[level]       = RA;
//    mlapiRP_[level]       = RP;
//    mlapiPR_[level]       = PR;
    mlapiRmod_[level]     = R;  
    mlapiPmod_[level]     = P;
    mlapiA_[level+1]      = C;
    mlapiS_[level]        = S;

    // break if coarsest level is below specified size
    if (C.GetNumGlobalRows() <= maxcoarsesize)
    {
      ++level;
      break;
    }
  
  } // for (level=1; level<maxlevels-1; ++level)
  
  // set coarse solver
  if (Comm().MyPID()==0)
  {
    ML_print_line("-", 78);
    cout << "VM3 solver : creating coarse solver level " << level << endl; 
    fflush(stdout);
  }
  S.Reshape(mlapiA_[level],coarsetype,mlparams_);

  mlapiS_[level] = S;
  
  // store number of levels
  maxlevels_ = level+1;

  iscomputed_ = true;
  return true;
}

