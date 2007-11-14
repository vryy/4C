#ifdef CCADISCRET

#include "vm3_solver.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 06/07|
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
 |  multigrid solver for VM3                                    vg 06/07|
 *----------------------------------------------------------------------*/
int VM3_Solver::Solve(const Epetra_Vector& B, Epetra_Vector& X, ParameterList& params)
{
    // do setup if not already done
  if (!iscomputed_) Compute();

  RCP<Epetra_Vector> x = LINALG::CreateVector(Acombined_->OperatorDomainMap(),true);
  RCP<Epetra_Vector> b = LINALG::CreateVector(Acombined_->OperatorRangeMap(),true);
  LINALG::Export(X,*x);
  LINALG::Export(B,*b);
  
  const Epetra_BlockMap& bmap = X.Map();
  Space space;
  space.Reshape(bmap.NumGlobalElements(),bmap.NumMyElements(),bmap.MyGlobalElements());
  
  MultiVector mvX(space,X.Pointers(),1);
  MultiVector mvB(space,B.Pointers(),1);
  
  MultiVector bcoarse;
  MultiVector xcoarse;
  bcoarse = Rtent_ * mvB;
  xcoarse = Rtent_ * mvX;
  
  RCP<Epetra_Vector> xcshifted = LINALG::CreateVector(*coarsermap_,true);
  RCP<Epetra_Vector> bcshifted = LINALG::CreateVector(*coarsermap_,true);
  const int mylength = xcshifted->MyLength();
  for (int i=0; i<mylength; ++i)
  {
    (*xcshifted)[i] = xcoarse(i,0);
    (*bcshifted)[i] = bcoarse(i,0);
  }
  LINALG::Export(*xcshifted,*x);
  LINALG::Export(*bcshifted,*b);
  
  RCP<ParameterList> rcpparams = rcp( new ParameterList(params));
  
  LINALG::Solver solver(rcpparams,Acombined_->RowMatrixRowMap().Comm(),NULL);
  solver.Solve(Acombined_,x,b,true,true);
  
  LINALG::Export(*x,*xcshifted);
  LINALG::Export(*x,X);
  for (int i=0; i<mylength; ++i)
    xcoarse(i,0) = (*xcshifted)[i];
 
  MultiVector x3h_h;
  x3h_h = Ptent_ * xcoarse;
  
  const int longlength = X.MyLength();
  for (int i=0; i<longlength; ++i)
    X[i] += x3h_h(i,0);

  return 0;

}

/*----------------------------------------------------------------------*
 |  compute the preconditioner (public)                         vg 06/07|
 *----------------------------------------------------------------------*/
bool VM3_Solver::Compute()
{
  // setup phase of multigrid
  iscomputed_ = false;

  // this is important to have!!!
  MLAPI::Init();

  // get parameters
  //int     maxlevels     = mlparams_.get("max levels",10);
  //int     maxcoarsesize = mlparams_.get("coarse: max size",10);
  double* nullspace     = mlparams_.get("null space: vectors",(double*)NULL);
  if (!nullspace) dserror("No nullspace supplied in parameter list");
  int     nsdim         = mlparams_.get("null space: dimension",1);
  //int     numpde        = mlparams_.get("PDE equations",1);
  //double  damping       = mlparams_.get("aggregation: damping factor",1.33);
  string  eigenanalysis = mlparams_.get("eigen-analysis: type", "Anorm");
  string  ptype         = mlparams_.get("prolongator: type","mod_full");
  string  coarsetype    = mlparams_.get("coarse: type","Amesos-KLU");
  string  smoothertype  = mlparams_.get("smoother: type","symmetric Gauss-Seidel");
  string  fsmoothertype = smoothertype;

  Space space(A_->RowMatrixRowMap());
  Operator mlapiA(space,space,A_.get(),false);
  Operator mlapiAplus(space,space,Aplus_.get(),false);

  // build nullspace;
  MultiVector NS;
  MultiVector NextNS;
  NS.Reshape(mlapiA.GetRangeSpace(),nsdim);
  if (nullspace)
  {
    const int length = NS.GetMyLength();
    for (int i=0; i<nsdim; ++i)
      for (int j=0; j<length; ++j)
        NS(j,i) = nullspace[i*length+j];
  }

  // get plain aggregation P and R
  Operator Ptent;
  Operator Rtent;
  GetPtent(mlapiA,mlparams_,NS,Ptent,NextNS);
  Rtent = GetTranspose(Ptent);

  // get coarse grid matrix K11 = R ( K+M ) P
  Operator K11;
  K11 = GetRAP(Rtent,mlapiA,Ptent);
  //cout << K11;

  Epetra_CrsMatrix* tmpcrs;
  int maxentries;
  double time;
  ML_Operator2EpetraCrsMatrix(K11.GetML_Operator(),tmpcrs,maxentries,true,time,0,true);
  RCP<Epetra_CrsMatrix> K11crs = rcp(tmpcrs);
  //cout << *K11crs;
  
  // get coarse grid matrix K12 = R ( K+M )
  Operator K12;
  K12 = Rtent * mlapiA;
  //cout << K12;
  ML_Operator2EpetraCrsMatrix(K12.GetML_Operator(),tmpcrs,maxentries,false,time,0,true);
  RCP<Epetra_CrsMatrix> K12crs = rcp(tmpcrs);
  //cout << *K12crs;

  // get fine grid matrix K21 = (K+M) * P
  Operator K21;
  K21 = mlapiA * Ptent;
  //cout << K21;
  ML_Operator2EpetraCrsMatrix(K21.GetML_Operator(),tmpcrs,maxentries,false,time,0,true);
  RCP<Epetra_CrsMatrix> K21crs = rcp(tmpcrs);
  //cout << *K21crs;
  
  // fine grid matrix is K22 = (K+M+M_fine);
  Operator K22 = mlapiAplus;
  
  // get row map of K22
  const Epetra_Map& k22rmap = Aplus_->RowMap();
  
  // get row map of K11
  const Epetra_Map& k11rmap = K11crs->RowMap();
  
  // build new map for row map of K11 that does not overlap with map of K22
  const int offset = k22rmap.MaxAllGID()+1;
  const int mygidsnewsize = k11rmap.NumMyElements();
  vector<int> mygidsnew(mygidsnewsize);
  for (int i=0; i<mygidsnewsize; ++i)
    mygidsnew[i] = k11rmap.GID(i) + offset;
  Epetra_Map k11rmapnew(-1,mygidsnewsize,&mygidsnew[0],0,k11rmap.Comm());

  // get combined row map of K11 and K22
  const int k22rmapsize = k22rmap.NumMyElements();
  const int mytotallength = mygidsnewsize + k22rmapsize;
  mygidsnew.resize(mytotallength);
  int last = 0;
  for (int i=0; i<k22rmapsize; ++i)
  {
    mygidsnew[i] = k22rmap.GID(i);
    last = i;
  }
  last +=1;
  const int k11rmapnewsize = k11rmapnew.NumMyElements();
  for (int i=0; i<k11rmapnewsize; ++i)
    mygidsnew[i+last] = k11rmapnew.GID(i);
  Epetra_Map kcombinedrmap(-1,mytotallength,&mygidsnew[0],0,k11rmap.Comm());
  
  // copy matrices K11crs and K12crs from row map k11rmap to k11rmapnew
  RCP<Epetra_CrsMatrix> K11new = LINALG::CreateMatrix(k11rmapnew,K11crs->MaxNumEntries()+10);
  RCP<Epetra_CrsMatrix> K12new = LINALG::CreateMatrix(k11rmapnew,K12crs->MaxNumEntries()+10);
  int length = max(K11crs->MaxNumEntries(),K12crs->MaxNumEntries());
  vector<int>    cindices(length);
  vector<double> values(length);
  const int k11rmapsize = k11rmap.NumMyElements();
  for (int i=0; i<k11rmapsize; ++i)
  {
    // new global row id
    const int grid = k11rmapnew.GID(i);
    int numindices = 0;
    int err = K11crs->ExtractMyRowCopy(i,length,numindices,&values[0],&cindices[0]);
    if (err<0) dserror("ExtractMyRowCopy returned %d",err);
    // new global column id gcid = (lcid -> gcid) + offset
    for (int j=0; j<numindices; ++j) cindices[j] = K11crs->ColMap().GID(cindices[j]) + offset;
    err = K11new->InsertGlobalValues(grid,numindices,&values[0],&cindices[0]);
    if (err<0) dserror("InsertGlobalValues returned %d",err);
    
    numindices = 0;
    err = K12crs->ExtractMyRowCopy(i,length,numindices,&values[0],&cindices[0]);
    // new global column id gcid = (lcid -> gcid)
    for (int j=0; j<numindices; ++j) cindices[j] = K12crs->ColMap().GID(cindices[j]);
    err = K12new->InsertGlobalValues(grid,numindices,&values[0],&cindices[0]);
    if (err<0) dserror("InsertGlobalValues returned %d",err);
  }
  K11new->FillComplete(k11rmapnew,k11rmapnew);
  K12new->FillComplete(k22rmap,k11rmapnew);

  // copy matrix K21 to new column map
  RCP<Epetra_CrsMatrix> K21new = LINALG::CreateMatrix(k22rmap,K21crs->MaxNumEntries()+10);
  length = K21crs->MaxNumEntries();
  cindices.resize(length);
  values.resize(length);
  for (int i=0; i<k22rmapsize; ++i)
  {
    int numindices=0;
    int err = K21crs->ExtractMyRowCopy(i,length,numindices,&values[0],&cindices[0]);
    int grid = k22rmap.GID(i);
    for (int j=0; j<numindices; ++j) cindices[j] = K21crs->ColMap().GID(cindices[j]) + offset;
    err = K21new->InsertGlobalValues(grid,numindices,&values[0],&cindices[0]);
    if (err<0) dserror("InsertGlobalValues returned %d",err);
  }
  K21new->FillComplete(k11rmapnew,k22rmap);

  // copy all matrices into one big matrix and store it
  int clength = K22.GetRowMatrix()->MaxNumEntries() + K11new->MaxNumEntries() + 100;
  RCP<Epetra_CrsMatrix> Kcombined = LINALG::CreateMatrix(kcombinedrmap,clength);
  LINALG::Add(*Aplus_,false,1.0,*Kcombined,0.0);
  LINALG::Add(*K21new,false,1.0,*Kcombined,1.0);
  LINALG::Add(*K12new,false,1.0,*Kcombined,1.0);
  LINALG::Add(*K11new,false,1.0,*Kcombined,1.0);
  Kcombined->FillComplete(kcombinedrmap,kcombinedrmap);
  Kcombined->OptimizeStorage();
  Acombined_ = Kcombined;

  // finally, we have to fix the nullspace in the parameter list to match the combined matrix dimension
  RCP<vector<double> > newnullsp = rcp(new vector<double>(nsdim*kcombinedrmap.NumMyElements()));
  for (int i=0; i<(int)newnullsp->size(); ++i) (*newnullsp)[i] = 1.0;
  mlparams_.set("null space: vectors",&(*newnullsp)[0]);
  mlparams_.set<RCP<vector<double> > >("nullspace",newnullsp);

  // store Ptent and Rtent
  Ptent_ = Ptent;
  Rtent_ = Rtent;
  
  // store new k11 row map
  coarsermap_ = rcp(new Epetra_Map(k11rmapnew));

  iscomputed_ = true;
  return true;
}



#endif
