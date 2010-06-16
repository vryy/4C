/*
 * amgpreconditioner.cpp
 *
 *  Created on: Jun 7, 2010
 *      Author: wiesner
 */

#ifdef CCADISCRET

#define WRITEOUTAGGREGATES

#include "Teuchos_TimeMonitor.hpp"

#include "amgpreconditioner.H"

#include "transfer_operator.H"
#include "transfer_operator_tentative.H"
#include "transfer_operator_saamg.H"
#include "transfer_operator_pgamg.H"

#include "aggregation_method_uncoupled.H"
#include "aggregation_method_ml.H"

#include "smoother.H"

LINALG::AMGPreconditioner::AMGPreconditioner(RCP<Epetra_Operator> A, const ParameterList& params, FILE* outfile)
:
  Epetra_Operator(),
  label_("LINALG::AMGPreconditioner"),
  params_(params),
  nlevel_(0),
  nmaxlevels_(0),
  nVerbose_(0)
{
  //SparseMatrix* tmp = dynamic_cast<SparseMatrix*>(A.get());
  SparseMatrix* tmp = new SparseMatrix((rcp_dynamic_cast<Epetra_CrsMatrix>(A)));
  if(!tmp) dserror("Expected Epetra_CrsMatrix");
  Ainput_ = rcp(tmp,false);
 // Ainput_ = rcp(new SparseMatrix(rcp_dynamic_cast<Epetra_CrsMatrix(A)));

  Setup();  // setup of multigrid hirarchy
}

void LINALG::AMGPreconditioner::Setup()
{
  // read parameters from parameter list
  int nmaxlevels_ = Params().get<int>("max levels",10);
  int maxcoarsesize = Params().get<int>("coarse: max size",1);
  int nsdim = Params().get<int>("null space: dimension",1);
  nVerbose_ = Params().get("output",0);

  params_.set("Unamalgamated BlockSize",nsdim);
  params_.set("PDE equations",nsdim);
  params_.set("ML output", nVerbose_);

  // prepare null space approximation
  const int length = Ainput_->RowMap().NumMyElements();
  const int nlnnode = length/nsdim;
  RCP<vector<double> > ns = rcp(new vector<double>(nsdim*length,0.0));
  for (int i=0; i<nlnnode; ++i)
  {
    (*ns)[i*nsdim] = 1.0;
    (*ns)[length+i*nsdim+1] = 1.0;
    (*ns)[2*length+i*nsdim+2] = 1.0;
    if(nsdim == 4)  // 3dim problems
      (*ns)[3*length+i*nsdim+3] = 1.0;    // note, this works only for constant number of dofs per node
  }


  params_.set("null space: vectors",&((*ns)[0])); // adapt default null space
  params_.remove("nullspace",false);
  RCP<Epetra_MultiVector> curNS = rcp(new Epetra_MultiVector(View,Ainput_->RowMap(),&((*ns)[0]),Ainput_->EpetraMatrix()->RowMatrixRowMap().NumMyElements(),nsdim));



  RCP<Epetra_MultiVector> nextNS = null;

  A_.resize(nmaxlevels_+1);    // transfer operators
  T_.resize(nmaxlevels_);    // transfer operators
  preS_.resize(nmaxlevels_);    // smoothers
  postS_.resize(nmaxlevels_);

  int curlevel = 0;
  A_[curlevel] = rcp(new SparseMatrix(*Ainput_,View));

#ifdef WRITEOUTAGGREGATES
    // plot out aggregates
    std::stringstream fileoutstream2;
    fileoutstream2 << "/home/wiesner/python/matrix" << curlevel << ".dat";
    LINALG::PrintMatrixInMatlabFormat(fileoutstream2.str(),*(A_[curlevel]->EpetraMatrix()));
#endif

  SmootherFactory smfac;

  for (curlevel = 0; curlevel < nmaxlevels_; ++curlevel)
  {
    /////////////////////////////////////////////////////////
    /////////////////////// AGGREGATION PROCESS
    RCP<Epetra_IntVector> aggs = null;

    ////////////// determine aggregates using the velocity block matrix A11_[curlevel]
    params_.set("phase 1: min nodes per aggregate", 9);
    params_.set("phase 1: max neighbour nodes", 2);
    params_.set("phase 2: node attachement scheme","MaxLink");
    params_.set("Current Level",curlevel);
    RCP<AggregationMethod> aggm = AggregationMethodFactory::Create("Uncoupled",NULL);
    int naggregates_local = 0;
    aggm->GetGlobalAggregates(A_[curlevel]->EpetraMatrix(),params_,aggs,naggregates_local,curNS);

#ifdef WRITEOUTAGGREGATES
    // plot out aggregates
    std::stringstream fileoutstream;
    fileoutstream << "/home/wiesner/python/aggregates" << curlevel << ".vel";
    aggm->PrintIntVectorInMatlabFormat(fileoutstream.str(),*aggs);
#endif

    /////////////////////////////////////////////////////////
    /////////////////////// CALCULATE TRANSFER OPERATORS

    ///////////// velocity transfer operators
    string ProlongSmoother = "PG-AMG"; //params_.get("prolongator smoother","PG-AMG"); TODO FIXME
    T_[curlevel] = TransferOperatorFactory::Create(ProlongSmoother,A_[curlevel],NULL);
    nextNS = T_[curlevel]->buildTransferOperators(aggs,naggregates_local,params_,curNS,0);


    A_[curlevel+1] = Multiply(T_[curlevel]->R(),*A_[curlevel],T_[curlevel]->P());

#ifdef WRITEOUTAGGREGATES
    // plot out aggregates
    std::stringstream fileoutstream2;
    fileoutstream2 << "/home/wiesner/python/matrix" << curlevel+1 << ".dat";
    LINALG::PrintMatrixInMatlabFormat(fileoutstream2.str(),*(A_[curlevel+1]->EpetraMatrix()));
#endif

    ParameterList IfpackParameters;
    IfpackParameters.set("relaxation: type","Gauss-Seidel");
    IfpackParameters.set("relaxation: sweeps",10);
    IfpackParameters.set("relaxation: damping factor", 0.66);
    IfpackParameters.set("relaxation: zero starting solution", false);
    IfpackParameters.set("fact: level-of-fill", 1);
    IfpackParameters.set("fact: ilut level-of-fill",0.0);
    IfpackParameters.set("amesos: solver type","Amesos_Umfpack");


    preS_[curlevel] = smfac.Create("ILU",A_[curlevel],IfpackParameters);
    postS_[curlevel] = preS_[curlevel];


    //////////////////// prepare variables for next aggregation level
    curNS = nextNS;
    nlevel_ = curlevel + 1;

    // abbruchbedingung
    if ((A_[curlevel+1]->EpetraMatrix()->NumGlobalRows() - aggm->getNumGlobalDirichletBlocks()*(nsdim)) < maxcoarsesize)
    {
      if(nVerbose_ > 5) cout << "dim A[" << curlevel+1 << "] < " << maxcoarsesize << ". -> end aggregation process" << endl;
      break;
    }

  }


  ParameterList IfpackParameters;
  IfpackParameters.set("amesos: solver type","Amesos_Umfpack");
  coarsestSmoother_ = smfac.Create("Amesos",A_[nlevel_],IfpackParameters);


}

int LINALG::AMGPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  RCP<Epetra_MultiVector> b = rcp(new Epetra_MultiVector(X));
  RCP<Epetra_MultiVector> x = rcp(new Epetra_MultiVector(Y));
  x->PutScalar(0.0);
  Vcycle(*b,*x,0);
  Y.Update(1.0,*x,0.0);
  return 0;

  // this is just an ILU with the same behaviour as the Vcycle -> very suspicious
  /*RCP<Epetra_MultiVector> b = rcp(new Epetra_MultiVector(X));
  RCP<Epetra_MultiVector> x = rcp(new Epetra_MultiVector(Y));
  x->PutScalar(0.0);
  preS_[0]->ApplyInverse(*b,*x);
  Y.Update(1.0,*x,0.0);
  return 0;*/

  // Use only the fine grid presmoother
  /*RCP<Epetra_MultiVector> bfine = rcp(new Epetra_MultiVector(X));
  RCP<Epetra_MultiVector> xfine = rcp(new Epetra_MultiVector(Y));

  preS_[0]->ApplyInverse(*bfine,*xfine);

  RCP<Epetra_MultiVector> rfine = rcp(new Epetra_MultiVector(X.Map(),1,true));
  A_[0]->Apply(*xfine,*rfine);
  rfine->Update(1.0,*bfine,-1.0);

  RCP<Epetra_MultiVector> bcoarse = rcp(new Epetra_MultiVector(T_[0]->Prolongator()->DomainMap(),1,true));
  T_[0]->Restrictor()->Apply(*rfine,*bcoarse);*/


}

#if 1
void LINALG::AMGPreconditioner::Vcycle(const Epetra_MultiVector& rhs, Epetra_MultiVector& sol, const int level) const
{
  // coarsest grid
  if(level == nlevel_)
  {
    coarsestSmoother_->ApplyInverse(rhs,sol);
    return;
  }

  ////////////// presmoothing
  preS_[level]->ApplyInverse(rhs,sol);

  ////////////// on finest level
  // calculate residual
  RCP<Epetra_MultiVector> res = rcp(new Epetra_MultiVector(rhs.Map(),1,true));
  if(level == 0)
  {
    A_[level]->Apply(sol,*res);
    res->Update(1.0,rhs,-1.0);
  }
  else
  {
    res->Update(1.0,rhs,0.0); // on carser grids, the rhs vector is the residual
  }

  /*double norm = 0.0;
  res->Norm2(&norm);
  cout << "PROC " << A_[level]->Comm().MyPID() << " Level " << level << " residual AFTER PRESMOOTHING norm " << norm << endl;*/

  ////////////// define vectors for coarse levels
  RCP<Epetra_MultiVector> rhs_c = rcp(new Epetra_MultiVector(T_[level]->OperatorDomainMap(),1,false));
  RCP<Epetra_MultiVector> sol_c = rcp(new Epetra_MultiVector(T_[level]->OperatorDomainMap(),1,true));

  ////////////// zero out coarse solution vector
  sol_c->PutScalar(0.0);

  ////////////// restrict residual to next coarser grid
  T_[level]->R().Apply(*res,*rhs_c);

  ////////////// recursive call of multigrid cycle
  Vcycle(*rhs_c,*sol_c,level + 1);

  ////////////// prolongate coarse solution
  T_[level]->P().Apply(*sol_c,*res);

  ////////////// update fine level solution
  sol.Update(1.0,*res,1.0);

  ////////////// postsmoothing
  if(level!=0)    // very suspicious, the postsmoother on level 0 somehow ignores the underlaying solution of the coarser grids.
    postS_[level]->ApplyInverse(rhs,sol);
  else
  {

    RCP<Epetra_MultiVector> rhs_tmp = rcp(new Epetra_MultiVector(rhs));
    RCP<Epetra_MultiVector> sol_tmp = rcp(new Epetra_MultiVector(sol)); // sol_tmp == sol!

    /* ONLY ILU */
    //A_[level]->Apply(sol,*rhs_tmp);
    //rhs_tmp->Update(-1.0,rhs,1.0); // rhs_tmp is difference of new rhs and old rhs


    // fine level postsmoother
    ParameterList IfpackParameters;
    IfpackParameters.set("relaxation: type","Gauss-Seidel");
    IfpackParameters.set("relaxation: sweeps",10);
    IfpackParameters.set("relaxation: damping factor", 0.66);
    IfpackParameters.set("relaxation: zero starting solution", false);
    IfpackParameters.set("fact: level-of-fill", 1);
    IfpackParameters.set("fact: ilut level-of-fill",0.0);
    IfpackParameters.set("amesos: solver type","Amesos_Umfpack");

    /*Smoother_Ifpack* prec = new Smoother_Ifpack("ILU", A_[0]->EpetraMatrix(),IfpackParameters);
    prec->ApplyInverse(*rhs_tmp,*sol_tmp);
    delete prec;*/

    SmootherFactory smfac;
    RCP<Smoother> prec = smfac.Create("ILU",A_[0],IfpackParameters);
    prec->ApplyInverse(*rhs_tmp,*sol_tmp);

    /*// this is working
    Ifpack factory;
    Ifpack_Preconditioner* prec = factory.Create("ILU",A_[0]->EpetraMatrix().get(),0);
    //Ifpack_Preconditioner* prec = factory.Create("point relaxation",A_[0]->EpetraMatrix().get(),0);

    prec->SetParameters(IfpackParameters);
    prec->Initialize();
    prec->Compute();
    if(prec->IsComputed()==false)
      dserror("prec not computed??");
    //prec->ApplyInverse(rhs,sol);
    prec->ApplyInverse(*rhs_tmp,*sol_tmp);  // sol tmp contains current solution (has no effect for ILU)

    sol.Update(-1.0,*sol_tmp,1.0);          // update solution sol ONLY ILU

    delete prec;*/

    /* my own small jacobi...
    for(int k=0; k<5;k++)
    {
      RCP<Epetra_MultiVector> tmpx = rcp(new Epetra_MultiVector(sol));
      A_[0]->Apply(sol,*tmpx);
      tmpx->Update(1.0,rhs,-1.0);
      sol.Update(0.66,*tmpx,1.0);
    }*/
  }

  /*res->PutScalar(0.0);
  A_[level]->Apply(sol,*res);
  res->Update(1.0,rhs,-1.0);
  res->Norm2(&norm);
  cout << "PROC " << A_[level]->Comm().MyPID() << " Level " << level << " residual AFTER POSTSMOOTHING norm " << norm << endl;
  */

  return;
}
#endif

#if 0
void LINALG::AMGPreconditioner::Vcycle(const Epetra_MultiVector& b_f, Epetra_MultiVector& x_f, const int level) const
{

  if(level == nlevel_/*-1*/)
  {
    // coarsest level
    coarsestSmoother_->ApplyInverse(b_f,x_f);

    return;
  }

  RCP<Epetra_MultiVector> r_c = rcp(new Epetra_MultiVector(T_[level]->OperatorDomainMap(),1,false));
  RCP<Epetra_MultiVector> x_c = rcp(new Epetra_MultiVector(T_[level]->OperatorDomainMap(),1,true));

  // presmoothing
  x_f.PutScalar(0.0);


  preS_[level]->ApplyInverse(b_f,x_f);


  // compute fine-level residual
  RCP<Epetra_MultiVector> r_f = rcp(new Epetra_MultiVector(b_f.Map(),1,true));
  A_[level]->Apply(x_f,*r_f);
  r_f->Update(1.0,b_f,-1.0);

  double norm = 0.0;
  r_f->Norm2(&norm);
  cout << "PROC " << A_[level]->Comm().MyPID() << " Level " << level << " residual AFTER PRESMOOTHING norm " << norm << endl;


  // restrict residual to next coarser grid
  T_[level]->R().Apply(*r_f,*r_c);

  // solve coarser problem
  Vcycle(*r_c,*x_c,level+1);

  // prolongate correction
  //r_f->PutScalar(0.0);
  RCP<Epetra_MultiVector> tmp = rcp(new Epetra_MultiVector(r_f->Map(),1,true));
  T_[level]->P().Apply(*x_c,*tmp);
  x_f.Update(1.0,*tmp,1.0);

  // postsmoothing
  postS_[level]->ApplyInverse(b_f,x_f);

  r_f->PutScalar(0.0);
  A_[level]->Apply(x_f,*r_f);
  r_f->Update(1.0,b_f,-1.0);
  r_f->Norm2(&norm);
  cout << "PROC " << A_[level]->Comm().MyPID() << " Level " << level << " residual AFTER POSTSMOOTHING norm " << norm << endl;

}
#endif

///////////////////////////////////////////////////////////////////
RCP<LINALG::SparseMatrix> LINALG::AMGPreconditioner::Multiply(const SparseMatrix& A, const SparseMatrix& B, const SparseMatrix& C, bool bComplete)
{
  TEUCHOS_FUNC_TIME_MONITOR("AMGPreconditioner::Multiply (with MLMultiply)");

  RCP<SparseMatrix> tmp = LINALG::MLMultiply(B,C,true);
  return LINALG::MLMultiply(A,*tmp,bComplete);
}



#endif
