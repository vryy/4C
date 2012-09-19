/*
 * amgpreconditioner.cpp
 *
 *  Created on: Jun 7, 2010
 *      Author: wiesner
 */



//#define WRITEOUTAGGREGATES
//#define ANALYSIS

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
  SparseMatrix* tmp = new SparseMatrix((rcp_dynamic_cast<Epetra_CrsMatrix>(A)));
  if(!tmp) dserror("Expected Epetra_CrsMatrix");
  Ainput_ = rcp(tmp,true);

  Setup();  // setup of multigrid hirarchy
}

void LINALG::AMGPreconditioner::Setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("AMGPreconditioner::Setup");

//  Epetra_Time time(A()->Comm());

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

  //cout << params_ << endl;

  for (curlevel = 0; curlevel < nmaxlevels_; ++curlevel)
  {
    /////////////////////////////////////////////////////////
    /////////////////////// AGGREGATION PROCESS
    RCP<Epetra_IntVector> aggs = null;

    ////////////// determine aggregates using the velocity block matrix A11_[curlevel]
    // set parameters for aggregation process
    params_.set("phase 1: min nodes per aggregate", params_.get("aggregation: nodes per aggregate",9));
    params_.set("phase 1: max neighbour nodes", 2);
    params_.set("phase 2: node attachement scheme","MaxLink");
    params_.set("Current Level",curlevel);
    if(params_.get("aggregation: threshold",0.0) != 0.0)
    {
      params_.set("aggregation method","anisotropic aggregation");
      params_.set("anisotropic aggregation: epsilon",params_.get("aggregation: threshold",0.3));
    }
    else
      params_.set("aggregation method","isotropic aggregation");

//    time.ResetStartTime();
    RCP<AggregationMethod> aggm = AggregationMethodFactory::Create("Uncoupled",NULL);
    int naggregates_local = 0;
    aggm->GetGlobalAggregates(A_[curlevel]->EpetraMatrix(),params_,aggs,naggregates_local,curNS);
//    cout << "Level " << curlevel << " AggregationTime: " << time.ElapsedTime() << endl; time.ResetStartTime();

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

//    cout << "Level " << curlevel << " PG-AMG Time " << time.ElapsedTime() << endl; time.ResetStartTime();

    A_[curlevel+1] = Multiply(T_[curlevel]->R(),*A_[curlevel],T_[curlevel]->P());

//    cout << "Level " << curlevel << " RAP time: " << time.ElapsedTime() << endl; time.ResetStartTime();

#ifdef WRITEOUTAGGREGATES
    // plot out aggregates
    std::stringstream fileoutstream2;
    fileoutstream2 << "/home/wiesner/python/matrix" << curlevel+1 << ".dat";
    LINALG::PrintMatrixInMatlabFormat(fileoutstream2.str(),*(A_[curlevel+1]->EpetraMatrix()));
#endif

    // extract smoother information
    std::stringstream strstream;
    strstream << "smoother: list (level " << curlevel << ")";
    ParameterList SmootherParams = params_.sublist(strstream.str());
    if(SmootherParams.get("smoother: type","Jacobi")=="Jacobi" ||
       SmootherParams.get("smoother: type","Jacobi")=="Gauss-Seidel")
    {
      ParameterList IfpackParameters;
      IfpackParameters.set("relaxation: type",SmootherParams.get("smoother: type","Jacobi"));
      IfpackParameters.set("relaxation: sweeps",SmootherParams.get("smoother: sweeps",1));
      IfpackParameters.set("relaxation: damping factor", SmootherParams.get("smoother: damping factor",0.66));
      IfpackParameters.set("relaxation: zero starting solution", false);
      IfpackParameters.set("schwarz: reordering type", "none");
      //IfpackParameters.set("partitioner: local parts",A_[curlevel]->EpetraMatrix()->NumMyRows()/3);

      preS_[curlevel] = smfac.Create("point relaxation"/*"block relaxation"*/,A_[curlevel],IfpackParameters);
      postS_[curlevel] = preS_[curlevel];

    }
    else if(SmootherParams.get("smoother: type","Jacobi") == "IFPACK")
    {
      if(SmootherParams.get("smoother: ifpack type","ILU") != "ILU") dserror("we expect an ILU smoother. others are not supported");

      ParameterList IfpackParameters;

      IfpackParameters.set("fact: level-of-fill",(int)SmootherParams.get("smoother: ifpack level-of-fill",1.0));
      IfpackParameters.set("fact: ilut level-of-fill",0.0);
      IfpackParameters.set("smoother: ifpack overlap", SmootherParams.get("smoother: ifpack overlap",0));
      IfpackParameters.set("schwarz: combine mode", params_.sublist("smoother: ifpack list").get("schwarz: combine mode","Zero"));
      IfpackParameters.set("schwarz: reordering type", params_.sublist("smoother: ifpack list").get("schwarz: reordering type","rcm"));
      IfpackParameters.set("partitioner: overlap", params_.sublist("smoother: ifpack list").get("partitioner: overlap",0));

      preS_[curlevel] = smfac.Create("ILU",A_[curlevel],IfpackParameters);
      postS_[curlevel] = preS_[curlevel];
    }
    else if(SmootherParams.get("smoother: type","Jacobi") == "Amesos-UMFPACK")
    {
      ParameterList IfpackParameters;
      IfpackParameters.set("amesos: solver type","Amesos_Umfpack");

      preS_[curlevel] = smfac.Create("Amesos",A_[curlevel],IfpackParameters);
      postS_[curlevel] = preS_[curlevel];
    }
    else
    {
      dserror("smoother type unknown");
    }

//    cout << "Level " << curlevel << " gen smoothers time: " << time.ElapsedTime() << endl;

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


  if(params_.get("coarse: type","Jacobi")=="Jacobi" ||
     params_.get("coarse: type","Jacobi")=="Gauss-Seidel")
  {
    ParameterList IfpackParameters;
    IfpackParameters.set("relaxation: type",params_.get("coarse: type","Jacobi"));
    IfpackParameters.set("relaxation: sweeps",params_.get("coarse: sweeps",1));
    IfpackParameters.set("relaxation: damping factor", params_.get("coarse: damping factor",0.66));
    IfpackParameters.set("schwarz: reordering type", "none");
    IfpackParameters.set("relaxation: zero starting solution", false);
    //IfpackParameters.set("partitioner: local parts",A_[nlevel_]->EpetraMatrix()->NumMyRows()/3);  // TODO fix me!


    coarsestSmoother_ = smfac.Create("point relaxation"/*"block relaxation"*/,A_[nlevel_],IfpackParameters);
  }
  else if(params_.get("coarse: type","Jacobi") == "IFPACK")
  {
    if(params_.get("coarse: ifpack type","ILU") != "ILU") dserror("we expect an ILU smoother. others are not supported");

    ParameterList IfpackParameters;
    IfpackParameters.set("fact: level-of-fill", params_.get("coarse: ifpack level-of-fill",1));
    IfpackParameters.set("fact: ilut level-of-fill",0.0);
    IfpackParameters.set("smoother: ifpack overlap", params_.get("coarse: ifpack overlap",0));
    IfpackParameters.set("schwarz: combine mode", params_.sublist("smoother: ifpack list").get("schwarz: combine mode","Add"));
    IfpackParameters.set("schwarz: reordering type", params_.sublist("smoother: ifpack list").get("schwarz: reordering type","rcm"));
    IfpackParameters.set("partitioner: overlap", params_.sublist("smoother: ifpack list").get("partitioner: overlap",0));

    coarsestSmoother_ = smfac.Create("ILU",A_[nlevel_],IfpackParameters);
  }
  else if(params_.get("coarse: type","Jacobi") == "Amesos-UMFPACK")
  {
    ParameterList IfpackParameters;
    IfpackParameters.set("amesos: solver type","Amesos_Umfpack");

    coarsestSmoother_ = smfac.Create("Amesos",A_[nlevel_],IfpackParameters);
  }
  else
  {
    dserror("smoother type unknown");
  }


}

int LINALG::AMGPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  TEUCHOS_FUNC_TIME_MONITOR("AMGPreconditioner::ApplyInverse");

  RCP<Epetra_MultiVector> b = rcp(new Epetra_MultiVector(X));
  RCP<Epetra_MultiVector> x = rcp(new Epetra_MultiVector(Y));
  x->PutScalar(0.0);
  Vcycle(*b,*x,0);
  Y.Update(1.0,*x,0.0);
  return 0;
}

void LINALG::AMGPreconditioner::Vcycle(const Epetra_MultiVector& rhs, Epetra_MultiVector& sol, const int level) const
{

  // coarsest grid
  if(level == nlevel_)
  {
    coarsestSmoother_->ApplyInverse(rhs,sol);
    return;
  }

#ifdef ANALYSIS
#ifdef DEBUG
  RCP<Epetra_Vector> resXX = rcp(new Epetra_Vector(sol.Map(),true));
  A_[level]->Apply(sol,*resXX);
  resXX->Update(1.0,rhs,-1.0);

  if(A_[level]->Comm().NumProc() == 1)
  {
    // runs only in serial
    FILE* fp = fopen("analysis.m","w");
    int length = A_[level]->RowMap().NumGlobalPoints();
    fprintf(fp,"A=sparse(%d,%d);\n",length,length);


    // loop over local rows of DinvAP0
    for(int row=0; row<A_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A_[level]->EpetraMatrix()->GRID(row);
      int nnz = A_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    fprintf(fp,"rhs = zeros(%d,1);\n",length);

    double * MyGlobalElements1 = resXX->Values();
    for(int i=0; i<resXX->Map().NumMyElements(); i++)
    {
      int grid = resXX->Map().GID(i);
      double vv = MyGlobalElements1[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,vv);
    }

    fclose(fp);

    cout << "** BEFORE presmoothing (level " << level << ") **" << endl << "Press y when you're done" << endl;
    char instring[100];
    scanf("%s",instring);
  }
#endif
#endif

  ////////////// presmoothing
  preS_[level]->ApplyInverse(rhs,sol);

  ////////////// on finest level
  // calculate residual
  RCP<Epetra_MultiVector> res = rcp(new Epetra_MultiVector(rhs.Map(),1,true));
  //if(level == 0)  // don't distinguish between levels! this is the best version...
  {
    A_[level]->Apply(sol,*res);
    res->Update(1.0,rhs,-1.0);
  }
  //else
  {
  //  res->Update(1.0,rhs,0.0); // on carser grids, the rhs vector is the residual
  }

#ifdef ANALYSIS
#ifdef DEBUG

  if(A_[level]->Comm().NumProc() == 1)
  {
    // runs only in serial
    FILE* fp = fopen("analysis.m","w");
    int length = A_[level]->RowMap().NumGlobalPoints();
    fprintf(fp,"A=sparse(%d,%d);\n",length,length);


    // loop over local rows of DinvAP0
    for(int row=0; row<A_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A_[level]->EpetraMatrix()->GRID(row);
      int nnz = A_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    fprintf(fp,"rhs = zeros(%d,1);\n",length);

    double * MyGlobalElements1 = res->Values();
    for(int i=0; i<res->Map().NumMyElements(); i++)
    {
      int grid = res->Map().GID(i);
      double vv = MyGlobalElements1[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,vv);
    }

    fclose(fp);

    cout << "** AFTER presmoothing (level " << level << ") **" << endl << "Press y when you're done" << endl;
    char instring[100];
    scanf("%s",instring);
  }
#endif
#endif

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

#ifdef ANALYSIS
#ifdef DEBUG
  resXX->PutScalar(0.0);
  A_[level]->Apply(sol,*resXX);
  resXX->Update(1.0,rhs,-1.0);

  if(A_[level]->Comm().NumProc() == 1)
  {
    // runs only in serial
    FILE* fp = fopen("analysis.m","w");
    int length = A_[level]->RowMap().NumGlobalPoints();
    fprintf(fp,"A=sparse(%d,%d);\n",length,length);


    // loop over local rows of DinvAP0
    for(int row=0; row<A_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A_[level]->EpetraMatrix()->GRID(row);
      int nnz = A_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    fprintf(fp,"rhs = zeros(%d,1);\n",length);

    double * MyGlobalElements1 = resXX->Values();
    for(int i=0; i<resXX->Map().NumMyElements(); i++)
    {
      int grid = resXX->Map().GID(i);
      double vv = MyGlobalElements1[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,vv);
    }

    fclose(fp);

    cout << "** AFTER coarse grid correction (level " << level << ") **" << endl << "Press y when you're done" << endl;
    char instring[100];
    scanf("%s",instring);
  }
#endif
#endif

  ////////////// postsmoothing
  postS_[level]->ApplyInverse(rhs,sol);

#ifdef ANALYSIS
#ifdef DEBUG
  resXX->PutScalar(0.0);
  A_[level]->Apply(sol,*resXX);
  resXX->Update(1.0,rhs,-1.0);

  if(A_[level]->Comm().NumProc() == 1)
  {
    // runs only in serial
    FILE* fp = fopen("analysis.m","w");
    int length = A_[level]->RowMap().NumGlobalPoints();
    fprintf(fp,"A=sparse(%d,%d);\n",length,length);


    // loop over local rows of DinvAP0
    for(int row=0; row<A_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A_[level]->EpetraMatrix()->GRID(row);
      int nnz = A_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    fprintf(fp,"rhs = zeros(%d,1);\n",length);

    double * MyGlobalElements1 = resXX->Values();
    for(int i=0; i<resXX->Map().NumMyElements(); i++)
    {
      int grid = resXX->Map().GID(i);
      double vv = MyGlobalElements1[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,vv);
    }

    fclose(fp);

    cout << "** AFTER postsmoothing (level " << level << ") **" << endl << "Press y when you're done" << endl;
    char instring[100];
    scanf("%s",instring);
  }
#endif
#endif
  return;
}

///////////////////////////////////////////////////////////////////
RCP<LINALG::SparseMatrix> LINALG::AMGPreconditioner::Multiply(const SparseMatrix& A, const SparseMatrix& B, const SparseMatrix& C, bool bComplete)
{
  TEUCHOS_FUNC_TIME_MONITOR("AMGPreconditioner::Multiply (with MLMultiply)");

  RCP<SparseMatrix> tmp = LINALG::MLMultiply(B,C,true);
  return LINALG::MLMultiply(A,*tmp,bComplete);
}



