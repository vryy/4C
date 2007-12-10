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
rplus_(rplus),
r_(r),
sol_(sol),
dbctoggle_(dbctoggle)
{
  if (compute_) Compute();
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 06/07|
 *----------------------------------------------------------------------*/
void VM3_Solver::Separate(RefCountPtr<Epetra_CrsMatrix>& Sep,
                     RefCountPtr<Epetra_CrsMatrix>& Aplus,
                     RefCountPtr<Epetra_CrsMatrix>& A,
                     RefCountPtr<Epetra_Vector>& rplus,
                     RefCountPtr<Epetra_Vector>& r,
                     RefCountPtr<Epetra_Vector>& sol,
                     bool increm)
{
  // pre- and post-multiply scale-separating operator matrix Sep
  Aplus_ = LINALG::MatMatMult(*Aplus_,false,*Sep_,false);
  Aplus_ = LINALG::MatMatMult(*Sep_,true,*Aplus_,false);

  // add the fine-scale matrix to obtain complete matrix
  LINALG::Add(*Aplus_,false,1.0,*A_,1.0);

  // compute the additional rhs and add it to existing rhs for incremental formul.
  if (increm)
  {
    Aplus_->Multiply(false,*sol_,*rplus_);
    r_->Update(-1.0,*rplus_,1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 06/07|
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
  //cout << *Sep;
  if (!Sep_->RowMap().SameAs(Aplus_->RowMap())) dserror("rowmap not equal");
  if (!Sep_->RangeMap().SameAs(Aplus_->RangeMap())) dserror("rangemap not equal");
  if (!Sep_->DomainMap().SameAs(Aplus_->DomainMap())) dserror("domainmap not equal");
  //if (!Sep->ColMap().SameAs(Aplus_->ColMap())) dserror("colmap not equal");

  return true;
}


#endif
