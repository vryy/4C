/*!----------------------------------------------------------------------
\file avms_solver.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "avms_solver.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_solver.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                               vg 06/07|
 *----------------------------------------------------------------------*/
AVMS_Solver::AVMS_Solver(RCP<Epetra_CrsMatrix> Aforfine,
                         RCP<Epetra_CrsMatrix> Aforcoarse,
                         const RCP<Epetra_Vector> dbctoggle,
                         ParameterList& mlparams) :
mlparams_(mlparams),
dbctoggle_(dbctoggle)
{
  Compute(Aforfine,Aforcoarse);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                               vg 06/07|
 *----------------------------------------------------------------------*/
AVMS_Solver::~AVMS_Solver()
{
  return;
}

/*----------------------------------------------------------------------*
 |  compute the preconditioner (private)                       gee 02/08|
 *----------------------------------------------------------------------*/
bool AVMS_Solver::Compute(RCP<Epetra_CrsMatrix> Aforfine,
                          RCP<Epetra_CrsMatrix> Aforcoarse)
{
  // this is important to have!!!
  MLAPI::Init();

  const Epetra_Vector& dbct = *dbctoggle_;

  // get nullspace parameters
  double* nullspace = mlparams_.get("null space: vectors",(double*)NULL);
  if (!nullspace) dserror("No nullspace supplied in parameter list");
  int nsdim = mlparams_.get("null space: dimension",1);

  // modify nullspace to ensure that DBC are fully taken into account
  if (nullspace)
  {
    const int length = Aforfine->RowMap().NumMyElements();
    for (int i=0; i<nsdim; ++i)
      for (int j=0; j<length; ++j)
        if (dbct[j]!=0.0) nullspace[i*length+j] = 0.0;
  }

  // get plain aggregation Ptent and Rtent
  RCP<Epetra_CrsMatrix>   Ptent;
  RCP<Epetra_MultiVector> nextNS;
  int offset = Aforcoarse->RangeMap().MaxAllGID() + 1;
  GetPtent(*Aforcoarse,mlparams_,nullspace,Ptent,nextNS,offset);

  // make K12 = P^T A
  RCP<Epetra_CrsMatrix> K12 = LINALG::Multiply(Ptent,true,Aforcoarse,false);

  // make K21 = A P
  RCP<Epetra_CrsMatrix> K21 = LINALG::Multiply(Aforcoarse,false,Ptent,false); 

  // get coarse grid matrix K11 = R ( K+M ) P = R K21
  RCP<Epetra_CrsMatrix> K11 = LINALG::Multiply(Ptent,true,K21,false);

  // rename Aplus_ into K22 for now
  RCP<Epetra_CrsMatrix> K22 = Aforfine;

  // merge maps of fine and coarse to build complete system
  int length = K22->RangeMap().NumMyElements() + K11->RangeMap().NumMyElements();
  vector<int> kgid(length);
  int count=0;
  for (int i=0; i<K22->RangeMap().NumMyElements(); ++i)
    kgid[count++] = K22->RangeMap().GID(i);
  for (int i=0; i<K11->RangeMap().NumMyElements(); ++i)
    kgid[count++] = K11->RangeMap().GID(i);
  if (count!=length) dserror("Dimension mismatch %d != %d",count,length);
  Epetra_Map kmap(-1,length,&kgid[0],0,K22->Comm());

  // create combined matrix and fill it
  Acombined_ = LINALG::CreateMatrix(kmap,K22->MaxNumEntries());
  LINALG::Add(*K22,false,1.0,*Acombined_,0.0);
  LINALG::Add(*K11,false,1.0,*Acombined_,1.0);
  LINALG::Add(*K12,false,1.0,*Acombined_,1.0);
  LINALG::Add(*K21,false,1.0,*Acombined_,1.0);
  int err = Acombined_->FillComplete(kmap,kmap);
  if (err) dserror("Epetra_CrsMatrix::FillComplete returned %d",err);
  err = Acombined_->OptimizeStorage();
  if (err) dserror("Epetra_CrsMatrix::OptimizaStorage returned %d",err);

  // finally, we have to fix the nullspace in the parameter list to 
  // match the combined matrix dimension
  if (nextNS->NumVectors() != nsdim) 
    dserror("Nullspace dimension mismatch between coarse and fine");
  RCP<vector<double> > newnullsp = rcp(new vector<double>(nsdim*kmap.NumMyElements()));
  count=0;
  for (int i=0; i<K22->RowMap().NumMyElements()*nsdim; ++i)
    (*newnullsp)[count++] = nullspace[i];
  for (int j=0; j<nsdim; ++j)
    for (int i=0; i<Ptent->DomainMap().NumMyElements(); ++i)
    {
      (*newnullsp)[count] = (*(*nextNS)(j))[i];
      ++count;
    }
  mlparams_.set("null space: vectors",&(*newnullsp)[0]);
  mlparams_.set<RCP<vector<double> > >("nullspace",newnullsp);

  // store Ptent
  P_ = Ptent;

  return true;
}

/*----------------------------------------------------------------------*
 |  multigrid solver for AVMS                                  gee 02/08|
 *----------------------------------------------------------------------*/
int AVMS_Solver::Solve(const Epetra_Vector& B, Epetra_Vector& X, 
                       ParameterList& params)
{
  // create combined vectors and export fine level guess and rhs to them
  RCP<Epetra_Vector> x = LINALG::CreateVector(Acombined_->DomainMap(),true);
  RCP<Epetra_Vector> b = LINALG::CreateVector(Acombined_->RangeMap(),true);
  LINALG::Export(X,*x);
  LINALG::Export(B,*b);

  // restrict initial guess and right hand side
  RCP<Epetra_Vector> xc = LINALG::CreateVector(P_->DomainMap(),false);
  RCP<Epetra_Vector> bc = LINALG::CreateVector(P_->DomainMap(),false);
  P_->Multiply(true,*x,*xc);
  P_->Multiply(true,*b,*bc);

  // export coarse values to combined vectors
  LINALG::Export(*xc,*x);
  LINALG::Export(*bc,*b);
  bc = null;
  
  // call solver
  {
    RCP<ParameterList> rcpparams = rcp( new ParameterList(params));
    LINALG::Solver solver(rcpparams,Acombined_->Comm(),NULL);
    solver.Solve(Acombined_,x,b,true,true);
    b = null;
  }
  
  // export coarse part of solution to coarse vector and fine part to fine vector
  LINALG::Export(*x,*xc);
  LINALG::Export(*x,X);
  x = null;
  
  // prolongate coarse part and add to fine part
  RCP<Epetra_Vector> xfine = LINALG::CreateVector(P_->RangeMap(),false);
  P_->Multiply(false,*xc,*xfine);
  xc = null;
  X.Update(1.0,*xfine,1.0);

  return 0;
}


#endif
