#ifdef CCADISCRET

#include "vm3_solver.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 06/07|
 *----------------------------------------------------------------------*/
VM3_Solver::VM3_Solver(RefCountPtr<Epetra_CrsMatrix>& Sep,
                       RefCountPtr<Epetra_CrsMatrix>& Aplus,
                       RefCountPtr<Epetra_CrsMatrix>& A,
                       RefCountPtr<Epetra_Vector>& sugrvisc,
                       RefCountPtr<Epetra_Vector>& rplus,
                       RefCountPtr<Epetra_Vector>& r,
                       RefCountPtr<Epetra_Vector>& sol,
                       const RefCountPtr<Epetra_Vector> dbctoggle,
                       ParameterList& mlparams,
                       bool compute) :
compute_(compute),
mlparams_(mlparams),
Sep_(Sep),
Aplus_(Aplus),
A_(A),
sugrvisc_(sugrvisc),
rplus_(rplus),
r_(r),
sol_(sol),
dbctoggle_(dbctoggle)
{
  if (compute_) Compute();
  return;
}

/*----------------------------------------------------------------------*
 |  scale-separating routine                                    vg 12/07|
 |  precomputation of unscaled S^T*M*S                                  |
 |  (called in the first timestep)                                      |
 *----------------------------------------------------------------------*/
void VM3_Solver::Separate(RefCountPtr<Epetra_CrsMatrix>& Sep,
                          RefCountPtr<Epetra_CrsMatrix>& Aplus)
{
  // pre- and post-multiply M by scale-separating operator matrix Sep
  Aplus_ = LINALG::MatMatMult(*Aplus_,false,*Sep_,false);
  Aplus_ = LINALG::MatMatMult(*Sep_,true,*Aplus_,false);

  return;
}

/*----------------------------------------------------------------------*
 |  scaling routine                                             vg 12/07|
 |  scale precomput. matrix product by subgrid-viscosity-scaling vector |
 |  (called in every timestep)                                          |
 *----------------------------------------------------------------------*/
void VM3_Solver::Scale(RefCountPtr<Epetra_CrsMatrix>& Aplus,
                       RefCountPtr<Epetra_CrsMatrix>& A,
                       RefCountPtr<Epetra_Vector>& sugrvisc,
                       RefCountPtr<Epetra_Vector>& rplus,
                       RefCountPtr<Epetra_Vector>& r,
                       RefCountPtr<Epetra_Vector>& sol,
                       bool increm)
{
  // compute the additional rhs and add it to existing rhs for incremental formul.
  if (increm)
  {
    // add the subgrid-viscosity-scaled fine-scale matrix to obtain complete matrix
    //AddVecScalMult(*Aplus_,*sol_,*sugrvisc_,*A_,*rplus_);
    VecScalMult(*Aplus_,*sol_,*sugrvisc_,*rplus_);
    r_->Update(-1.0,*rplus_,1.0);
  }
  else
    // add the subgrid-viscosity-scaled fine-scale matrix to obtain complete matrix
    AddVecScal(*Aplus_,*sugrvisc_,*A_);


  return;
}


/*----------------------------------------------------------------------*
 |  Add a sparse matrix to another, the first one vector-scaled vg 12/07|
 |  B = B + A*V                                                         |
 *----------------------------------------------------------------------*/
void VM3_Solver::AddVecScal(const Epetra_CrsMatrix& A,
                            const Epetra_Vector& V,
                                  Epetra_CrsMatrix& B)
{
  if (!A.Filled()) dserror("FillComplete was not called on A");
  if (B.Filled()) dserror("FillComplete was called on B before");

  Epetra_CrsMatrix*               Aprime = NULL;
  Aprime = const_cast<Epetra_CrsMatrix*>(&A);

  //Loop over Aprime's rows and sum into
  int MaxNumEntries = EPETRA_MAX( Aprime->MaxNumEntries(), B.MaxNumEntries() );
  int NumEntries;
  vector<int>    Indices(MaxNumEntries);
  vector<double> Values(MaxNumEntries);

  const int NumMyRows = Aprime->NumMyRows();
  int Row, err, Col;
  for( int i = 0; i < NumMyRows; ++i )
  {
    Row = Aprime->GRID(i);
    int ierr = Aprime->ExtractGlobalRowCopy(Row,MaxNumEntries,NumEntries,&Values[0],&Indices[0]);
    if (ierr) dserror("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d",ierr);
    for( int j=0; j<NumEntries; ++j ) 
    {
      Col = Indices[j];
      Values[j] *= (sqrt(V[Row])*sqrt(V[Col]));
    }
    for (int j=0; j<NumEntries; ++j)
    {
      err = B.SumIntoGlobalValues(Row,1,&Values[j],&Indices[j]);
      if (err<0 || err==2)
        err = B.InsertGlobalValues(Row,1,&Values[j],&Indices[j]);
      if (err < 0)
        dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Add a sparse matrix to another, the first one vector-scaled vg 12/07|
 |  B = B + A*V                                                         |
 |  and calculate right hand side by multiplying vector-scaled (A*V)    |
 |  with solution vector S: R = (A*V)*S                                 |
 *----------------------------------------------------------------------*/
void VM3_Solver::AddVecScalMult(const Epetra_CrsMatrix& A,
                                const Epetra_Vector& S,
                                const Epetra_Vector& V,
                                      Epetra_CrsMatrix& B,
                                      Epetra_Vector& R)
{
  if (!A.Filled()) dserror("FillComplete was not called on A");
  if (B.Filled()) dserror("FillComplete was called on B before");

  Epetra_CrsMatrix*               Aprime = NULL;
  Aprime = const_cast<Epetra_CrsMatrix*>(&A);

  //Loop over Aprime's rows and sum into
  int MaxNumEntries = EPETRA_MAX( Aprime->MaxNumEntries(), B.MaxNumEntries() );
  int NumEntries;
  vector<int>    Indices(MaxNumEntries);
  vector<double> Values(MaxNumEntries);

  const int NumMyRows = Aprime->NumMyRows();
  int Row, err, Col;
  for( int i = 0; i < NumMyRows; ++i )
  {
    Row = Aprime->GRID(i);
    int ierr = Aprime->ExtractGlobalRowCopy(Row,MaxNumEntries,NumEntries,&Values[0],&Indices[0]);
    if (ierr) dserror("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d",ierr);
    for( int j=0; j<NumEntries; ++j ) 
    {
      Col = Indices[j];
      Values[j] *= (sqrt(V[Row])*sqrt(V[Col]));
      R[i] += Values[j]*S[j];
    }
    for (int j=0; j<NumEntries; ++j)
    {
      err = B.SumIntoGlobalValues(Row,1,&Values[j],&Indices[j]);
      if (err<0 || err==2)
        err = B.InsertGlobalValues(Row,1,&Values[j],&Indices[j]);
      if (err < 0)
        dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Calculate right hand side by multiplying vector-scaled (A*V)        |
 |  with solution vector S: R = (A*V)*S                                 |
 *----------------------------------------------------------------------*/
void VM3_Solver::VecScalMult(const Epetra_CrsMatrix& A,
                             const Epetra_Vector& S,
                             const Epetra_Vector& V,
                                   Epetra_Vector& R)
{
  Epetra_CrsMatrix*               Aprime = NULL;
  Aprime = const_cast<Epetra_CrsMatrix*>(&A);

  //Loop over Aprime's rows and sum into
  int MaxNumEntries = Aprime->MaxNumEntries();
  int NumEntries;
  vector<int>    Indices(MaxNumEntries);
  vector<double> Values(MaxNumEntries);

  const int NumMyRows = Aprime->NumMyRows();
  int Row, Col;
  for( int i = 0; i < NumMyRows; ++i )
  {
    Row = Aprime->GRID(i);
    int ierr = Aprime->ExtractGlobalRowCopy(Row,MaxNumEntries,NumEntries,&Values[0],&Indices[0]);
    if (ierr) dserror("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d",ierr);
    for( int j=0; j<NumEntries; ++j ) 
    {
      Col = Indices[j];
      Values[j] *= (sqrt(V[Row])*sqrt(V[Col]));
      R[i] += Values[j]*S[j];
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  routine for generating the scale-separation matrix S        vg 12/07|
 |  (will only be called in the first timestep)                         |
 *----------------------------------------------------------------------*/
bool VM3_Solver::Compute()
{
  // this is important to have!!!
  MLAPI::Init();

  const Epetra_Vector& dbct = *dbctoggle_;

  // get nullspace parameters
  double* nullspace     = mlparams_.get("null space: vectors",(double*)NULL);
  if (!nullspace) dserror("No nullspace supplied in parameter list");
  int     nsdim         = mlparams_.get("null space: dimension",1);

  Space space(Aplus_->RowMap());
  Operator mlapiAplus(space,space,Aplus_.get(),false);

  // build nullspace;
  MultiVector NS;
  MultiVector NextNS;
  NS.Reshape(mlapiAplus.GetRangeSpace(),nsdim);
  if (nullspace)
  {
    const int length = NS.GetMyLength();
    for (int i=0; i<nsdim; ++i)
     for (int j=0; j<length; ++j)
        if (dbct[j]!=1.0) NS(j,i) = nullspace[i*length+j];
  }

  // get plain aggregation P and R
  Operator Ptent;
  Operator Rtent;
  GetPtent(mlapiAplus,mlparams_,NS,Ptent,NextNS);
  Rtent = GetTranspose(Ptent);

  // get scale-separating operator S
  Operator S;
  Operator I = GetIdentity(mlapiAplus.GetDomainSpace(),mlapiAplus.GetRangeSpace());
  S = I - (Ptent * Rtent);

  // transfer scale-separating operator S from MLAPI to EpetraCrsMatrix
  Epetra_CrsMatrix* tmpcrs;
  int maxentries;
  double time;
  ML_Operator2EpetraCrsMatrix(S.GetML_Operator(),tmpcrs,maxentries,true,time,0,true);
  RCP<Epetra_CrsMatrix> tmprcp = rcp(tmpcrs);
  Epetra_Export exporter(tmprcp->RowMap(),Aplus_->RowMap());
  int err = Sep_->Export(*tmprcp,exporter,Insert);
  if (err) dserror("Export returned %d",err);
  if (!Sep_->Filled()) Sep_->FillComplete(Aplus_->DomainMap(),Aplus_->RangeMap());
  tmprcp = null;

  if (!Sep_->RowMap().SameAs(Aplus_->RowMap())) dserror("rowmap not equal");
  if (!Sep_->RangeMap().SameAs(Aplus_->RangeMap())) dserror("rangemap not equal");
  if (!Sep_->DomainMap().SameAs(Aplus_->DomainMap())) dserror("domainmap not equal");
  //if (!Sep->ColMap().SameAs(Aplus_->ColMap())) dserror("colmap not equal");

  return true;
}


#endif
