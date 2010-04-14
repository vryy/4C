/*
 * saddlepointpreconditioner.cpp
 *
 *  Created on: Feb 16, 2010
 *      Author: wiesner
 */

#ifdef CCADISCRET


#define WRITEOUTSTATISTICS
#undef WRITEOUTSYMMETRY // write out || A - A^T ||_F
#define WRITEOUTAGGREGATES

#include "linalg_sparsematrix.H"

#include "saddlepointpreconditioner.H"
#include "braesssarazin_smoother.H"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "MLAPI_Aggregation.h"

// includes for MLAPI functions
#include "ml_common.h"
#include "ml_include.h"
#include "ml_aggregate.h"
#include "ml_agg_METIS.h"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Error.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI.h"

#include "float.h" // for DBL_MAX and DBL_MIN

LINALG::SaddlePointPreconditioner::SaddlePointPreconditioner(RCP<Epetra_Operator> A, const ParameterList& params, const ParameterList& pressurelist, FILE* outfile)
: params_(params),
pressureparams_(pressurelist),
outfile_(outfile),
nVerbose_(0)
{
  Setup(A,params,pressurelist);
}

LINALG::SaddlePointPreconditioner::~SaddlePointPreconditioner()
{

}


int LINALG::SaddlePointPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // VCycle

  // note: Aztec might pass X and Y as physically identical objects,
  // so we better deep copy here

  RCP<LINALG::ANA::Vector> Xv = rcp(new LINALG::ANA::Vector(*mmex_.Map(0),false));
  RCP<LINALG::ANA::Vector> Xp = rcp(new LINALG::ANA::Vector(*mmex_.Map(1),false));

  RCP<LINALG::ANA::Vector> Yv = rcp(new LINALG::ANA::Vector(*mmex_.Map(0),false));
  RCP<LINALG::ANA::Vector> Yp = rcp(new LINALG::ANA::Vector(*mmex_.Map(1),false));


  // split vector using mmex_
  mmex_.ExtractVector(X,0,*Xv);
  mmex_.ExtractVector(X,1,*Xp);

  VCycle(*Xv,*Xp,*Yv,*Yp,0);

  mmex_.InsertVector(*Yv,0,Y);
  mmex_.InsertVector(*Yp,1,Y);

  return 0;
}

int LINALG::SaddlePointPreconditioner::VCycle(const Epetra_MultiVector& Xvel, const Epetra_MultiVector& Xpre, Epetra_MultiVector& Yvel, Epetra_MultiVector& Ypre, const int level) const
{
  // Y = A_^{-1} * X => solve A*Y = X

  if (level == nlevels_)
  {
    // coarsest level
    coarsestSmoother_->ApplyInverse(Xvel,Xpre,Yvel,Ypre);

    return 0;
  }

  // vectors for presmoothed solution
  RCP<Epetra_MultiVector> Zvel = rcp(new Epetra_MultiVector(Yvel.Map(),1,true));
  RCP<Epetra_MultiVector> Zpre = rcp(new Epetra_MultiVector(Ypre.Map(),1,true));

  // presmoothing
  if(bPresmoothing_) preS_[level]->ApplyInverse(Xvel,Xpre,*Zvel,*Zpre);  // rhs X is fix, initial solution Z = 0 (per definition, see above)
                                                      // note: ApplyInverse expects the "solution" and no solution increment "Delta Z"

  // calculate residual (fine grid)
  RCP<Epetra_Vector> velres = rcp(new Epetra_Vector(Yvel.Map(),true));
  RCP<Epetra_Vector> preres = rcp(new Epetra_Vector(Ypre.Map(),true));

  RCP<Epetra_Vector> vtemp = rcp(new Epetra_Vector(Yvel.Map(),true));
  RCP<Epetra_Vector> ptemp = rcp(new Epetra_Vector(Ypre.Map(),true));

  A11_[level]->Apply(*Zvel,*vtemp);
  A12_[level]->Apply(*Zpre,*velres);
  velres->Update(1.0,*vtemp,1.0);  // velres = + F Zvel + G Zpre
  velres->Update(1.0,Xvel,-1.0); // velres = Xvel - F Zvel - G Zpre

  A21_[level]->Apply(*Zvel,*ptemp);
  A22_[level]->Apply(*Zpre,*preres);
  preres->Update(1.0,*ptemp,1.0); // preres = + D Zvel + Z Zpre
  preres->Update(1.0,Xpre,-1.0); // preres = Xpre - D Zvel - Z Zpre

  // calculate coarse residual
  RCP<Epetra_Vector> velres_coarse = rcp(new Epetra_Vector(Rvel_[level]->RowMap(),true));
  RCP<Epetra_Vector> preres_coarse = rcp(new Epetra_Vector(Rpre_[level]->RowMap(),true));
  Rvel_[level]->Apply(*velres,*velres_coarse);
  Rpre_[level]->Apply(*preres,*preres_coarse);

  // define vector for coarse level solution
  RCP<Epetra_Vector> velsol_coarse = rcp(new Epetra_Vector(A11_[level+1]->RowMap(),true));
  RCP<Epetra_Vector> presol_coarse = rcp(new Epetra_Vector(A22_[level+1]->RowMap(),true));

  // call Vcycle recursively
  VCycle(*velres_coarse,*preres_coarse,*velsol_coarse,*presol_coarse,level+1);

  // define vectors for prolongated solution
  RCP<Epetra_Vector> velsol_prolongated = rcp(new Epetra_Vector(A11_[level]->RowMap(),true));
  RCP<Epetra_Vector> presol_prolongated = rcp(new Epetra_Vector(A22_[level]->RowMap(),true));

  // prolongate solution
  Pvel_[level]->Apply(*velsol_coarse,*velsol_prolongated);
  Ppre_[level]->Apply(*presol_coarse,*presol_prolongated);

  // update solution Zvel and Zpre for postsmoother
  Zvel->Update(1.0,*velsol_prolongated,1.0);
  Zpre->Update(1.0,*presol_prolongated,1.0);

  // postsmoothing
  if (bPostsmoothing_) postS_[level]->ApplyInverse(Xvel,Xpre,*Zvel,*Zpre); // rhs the same as for presmoothing, but better initial solution (Z)

  // write out solution
  Yvel.Update(1.0,*Zvel,0.0);
  Ypre.Update(1.0,*Zpre,0.0);

  return 0;
}

void LINALG::SaddlePointPreconditioner::Setup(RCP<Epetra_Operator> A,const ParameterList& origlist,const ParameterList& origplist)
{

#ifdef WRITEOUTSTATISTICS
  Epetra_Time ttt(A->Comm());
  ttt.ResetStartTime();
#endif




  // SETUP with SparseMatrix base class
  //////////////////// define some variables
  //const int myrank = A->Comm().MyPID();
  Epetra_Time time(A->Comm());
  const Epetra_Map& fullmap = A->OperatorRangeMap();
  const int         length  = fullmap.NumMyElements();
  nVerbose_ = 0;      // level of verbosity
  int ndofpernode = 0;// number of dofs per node
  int nv = 0;         // number of velocity dofs per node
  int np = 0;         // number of pressure dofs per node (1)
  int nlnode;         // number of nodes (local)


  Teuchos::RCP<Epetra_MultiVector> curvelNS = null;   // variables for null space
  Teuchos::RCP<Epetra_MultiVector> nextvelNS = null;
  Teuchos::RCP<Epetra_MultiVector> curpreNS = null;
  Teuchos::RCP<Epetra_MultiVector> nextpreNS = null;

  ///////////////// set parameter list
  RCP<ParameterList> spparams = rcp(new ParameterList());     // all paramaters
  RCP<ParameterList> velparams = rcp(new ParameterList());    // parameters (velocity specific)
  RCP<ParameterList> preparams = rcp(new ParameterList());    // parameters (pressure specific)

  // obtain common ML parameters from FLUID SOLVER block for coarsening from the dat file
  // we need at least "ML Parameters"."PDE equations" and "nullspace" information
  spparams->sublist("AMGBS Parameters") = params_.sublist("AMGBS Parameters"); // copy common parameters
  spparams->sublist("AMGBS Parameters").set("PDE equations",params_.sublist("ML Parameters").get("PDE equations",3));
  spparams->sublist("AMGBS Parameters").set("null space: add default vectors",params_.sublist("ML Parameters").get("null space: add default vectors",false));
  spparams->sublist("AMGBS Parameters").set("null space: dimension",params_.sublist("ML Parameters").get("null space: dimension",3));
  //spparams->sublist("AMGBS Parameters").set("null space: vectors",params_.sublist("ML Parameters").get("null space: vectors",NULL));
  spparams->sublist("AMGBS Parameters").set("ML output",spparams->sublist("AMGBS Parameters").get("output",0)); // set ML output
  spparams->sublist("AMGBS Parameters").remove("output");
  spparams->sublist("AMGBS Parameters").remove("smoother: type");  // we're using Braess-Sarazin only

  params_.remove("ML Parameters",false);  // now we don't need the ML Parameters any more


  /////////////////// prepare variables
  nmaxlevels_ = spparams->sublist("AMGBS Parameters").get("max levels",6) - 1;
  nlevels_ = 0;       // no levels defined
  bPresmoothing_ = false;   // get flags for pre- and postsmoothing
  bPostsmoothing_ = false;
  if(spparams->sublist("AMGBS Parameters").get("amgbs: smoother: pre or post","both") == "both" ||
     spparams->sublist("AMGBS Parameters").get("amgbs: smoother: pre or post","both") == "pre")
    bPresmoothing_ = true;
  if(spparams->sublist("AMGBS Parameters").get("amgbs: smoother: pre or post","both") == "both" ||
     spparams->sublist("AMGBS Parameters").get("amgbs: smoother: pre or post","both") == "post")
    bPostsmoothing_ = true;
  A11_.resize(nmaxlevels_+1);
  A12_.resize(nmaxlevels_+1);
  A21_.resize(nmaxlevels_+1);
  A22_.resize(nmaxlevels_+1);
  Pvel_.resize(nmaxlevels_);
  Ppre_.resize(nmaxlevels_);
  Rvel_.resize(nmaxlevels_);
  Rpre_.resize(nmaxlevels_);
  preS_.resize(nmaxlevels_);
  postS_.resize(nmaxlevels_);

  int nmaxcoarsedim = spparams->sublist("AMGBS Parameters").get("max coarse dimension",20);
  nVerbose_ = spparams->sublist("AMGBS Parameters").get("ML output",0);
  ndofpernode = spparams->sublist("AMGBS Parameters").get<int>("PDE equations",0);
  if(ndofpernode == 0) dserror("dof per node is zero -> error");

  nv       = ndofpernode-1;
  np       = 1;
  nlnode   = length / ndofpernode;



  /////////////////// transform Input matrix
  Ainput_ = rcp_dynamic_cast<BlockSparseMatrixBase>(A);
  if(Ainput_ != null)
  {
    mmex_ = Ainput_->RangeExtractor();
  }
  else
  {
    // get # dofs per node from params_ list and split row map
    time.ResetStartTime();
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
    //if (!myrank /*&& SIMPLER_TIMING*/) printf("--- Time to split map       %10.3E\n",time.ElapsedTime());
    time.ResetStartTime();
    // wrap matrix in SparseMatrix and split it into 2x2 BlockMatrix
    {
      SparseMatrix fullmatrix(rcp_dynamic_cast<Epetra_CrsMatrix>(A));
      Ainput_ = fullmatrix.Split<LINALG::DefaultBlockMatrixStrategy>(mmex_,mmex_);
      //if (!myrank /*&& SIMPLER_TIMING*/) printf("--- Time to split matrix    %10.3E\n",time.ElapsedTime());
      time.ResetStartTime();
      Ainput_->Complete();
      //if (!myrank /*&& SIMPLER_TIMING*/) printf("--- Time to complete matrix %10.3E\n",time.ElapsedTime());
      time.ResetStartTime();
    }
  }

  /////////////////// prepare null space for finest level (splitted into velocity and pressure part)

  // velocity part: fill in parameter list
  velparams->sublist("AMGBS Parameters") = spparams->sublist("AMGBS Parameters"); // copy common parameters
  velparams->sublist("AMGBS Parameters").set("PDE equations",nv);             // adapt nPDE (only velocity dofs)
  velparams->sublist("AMGBS Parameters").set("null space: dimension",nv);
  const int vlength = (*Ainput_)(0,0).RowMap().NumMyElements();
  RCP<vector<double> > vnewns = rcp(new vector<double>(nv*vlength,0.0));
  for (int i=0; i<nlnode; ++i)
  {
    (*vnewns)[i*nv] = 1.0;
    (*vnewns)[vlength+i*nv+1] = 1.0;
    if (nv>2) (*vnewns)[2*vlength+i*nv+2] = 1.0;
  }
  velparams->sublist("AMGBS Parameters").set("null space: vectors",&((*vnewns)[0])); // adapt default null space
  velparams->sublist("AMGBS Parameters").remove("nullspace",false);

  curvelNS = rcp(new Epetra_MultiVector(View,(*Ainput_)(0,0).RowMap(),&((*vnewns)[0]),(*Ainput_)(0,0).EpetraMatrix()->RowMatrixRowMap().NumMyElements(),nv));


  // pressure part: fill parameter list
  preparams->sublist("AMGBS Parameters") = spparams->sublist("AMGBS Parameters");
  preparams->sublist("AMGBS Parameters").set("PDE equations",1);               // adapt nPDE (only one pressure dof)
  preparams->sublist("AMGBS Parameters").set("null space: dimension", 1);
  const int plength = (*Ainput_)(1,1).RowMap().NumMyElements();
  RCP<vector<double> > pnewns = rcp(new vector<double>(plength,1.0));
  preparams->sublist("AMGBS Parameters").set("null space: vectors",&((*pnewns)[0]));
  preparams->sublist("AMGBS Parameters").remove("nullspace",false);

  curpreNS = rcp(new Epetra_MultiVector(View,(*Ainput_)(1,1).RowMap(),&((*pnewns)[0]),(*Ainput_)(1,1).EpetraMatrix()->RowMatrixRowMap().NumMyElements(),1));

  ////////////////// store level 0 matrices (finest level)
  int curlevel = 0;

  A11_[curlevel] = rcp(new SparseMatrix(Ainput_->Matrix(0,0),Copy));    // check me: copy or view only??
  A12_[curlevel] = rcp(new SparseMatrix(Ainput_->Matrix(0,1),Copy));
  A21_[curlevel] = rcp(new SparseMatrix(Ainput_->Matrix(1,0),Copy));
  A22_[curlevel] = rcp(new SparseMatrix(Ainput_->Matrix(1,1),Copy));

  MLAPI::Init();

  for (curlevel = 0; curlevel < nmaxlevels_; ++curlevel)
  {
    /////////////////////////////////////////////////////////
    /////////////////////// CALCULATE PTENT
    RCP<Epetra_IntVector> velaggs = rcp(new Epetra_IntVector(A11_[curlevel]->RowMap(),true));
    RCP<Epetra_IntVector> preaggs = rcp(new Epetra_IntVector(A22_[curlevel]->RowMap(),true));

    RCP<SparseMatrix> vel_Ptent   = null;  // these variables are filled by GetPtent
    RCP<SparseMatrix> pre_Ptent   = null;

    // determine aggregates using the velocity block matrix A11_[curlevel]
    int naggregates_local = 0;
    int naggregates = GetGlobalAggregates(*A11_[curlevel],velparams->sublist("AMGBS Parameters"),*curvelNS,*velaggs,naggregates_local);

    // transform tmp_velPtent to MLAPI Operator
    GetPtent(A11_[curlevel]->RowMap(),*velaggs,naggregates_local,velparams->sublist("AMGBS Parameters"),*curvelNS,vel_Ptent,nextvelNS,0);

    // transform vector with velocity aggregates to pressure block
    for(int i=0; i < preaggs->MyLength(); i++)
    {
      (*preaggs)[i] = (*velaggs)[i*nv];
    }

    // calculating Ptent for pressure block

    GetPtent(A22_[curlevel]->RowMap(),*preaggs,naggregates_local,preparams->sublist("AMGBS Parameters"),*curpreNS,pre_Ptent,nextpreNS,naggregates*nv);

#ifdef WRITEOUTAGGREGATES
//    std::ofstream fileout;
//    std::stringstream fileoutstream;
//    fileoutstream << "/home/wiesner/fluid/cubits/dc/aggregates/dc/aggregates" << curlevel << ".vel";
//    fileout.open(fileoutstream.str().c_str(),ios_base::out);
//    velaggs->Print(fileout);
//    fileout.flush();
//    fileout.close();
//
//    std::stringstream fileoutstreamp;
//    fileoutstreamp << "/home/wiesner/fluid/cubits/dc/aggregates/dc/aggregates" << curlevel << ".pre";
//    fileout.open(fileoutstreamp.str().c_str(),ios_base::out);
//    preaggs->Print(fileout);
//    fileout.flush();
//    fileout.close();


 /*   std::ofstream fileout2;
    std::stringstream fileoutstream2;
    fileoutstream2 << "/home/wiesner/Amat" << curlevel << ".txt";
    fileout2.open(fileoutstream2.str().c_str(),ios_base::out);
    fileout2 << *A11_[curlevel] << endl;
    fileout2.flush();
    fileout2.close();*/
#endif

    ////////////////////////////////////////////////
    /////////////////// CALCULATE RTENT

    // just transpose tentative prolongators (unsmoothed prolongators!!!)
    RCP<SparseMatrix> vel_Rtent = vel_Ptent->Transpose();
    RCP<SparseMatrix> pre_Rtent = pre_Ptent->Transpose();


    /////////////////////////// prolongator and restrictor smoothing

    RCP<SparseMatrix> vel_Psmoothed = null;
    RCP<SparseMatrix> pre_Psmoothed = null;
    RCP<SparseMatrix> vel_Rsmoothed = null;
    RCP<SparseMatrix> pre_Rsmoothed = null;

    // smooth velocity part
    string velProlongSmoother = velparams->sublist("AMGBS Parameters").get("amgbs: prolongator smoother (vel)","PA-AMG");
    if(velProlongSmoother == "SA-AMG") // SMOOTHED AGGREGATION
    {
      SA_AMG(A11_[curlevel],vel_Ptent,vel_Rtent,vel_Psmoothed,vel_Rsmoothed);

      Pvel_[curlevel] = vel_Psmoothed;
      Rvel_[curlevel] = vel_Rsmoothed;
    }
    else if(velProlongSmoother == "PG-AMG") // PETROV GALERKIN
    {
      PG_AMG(A11_[curlevel],vel_Ptent,vel_Rtent,vel_Psmoothed,vel_Rsmoothed);

      Pvel_[curlevel] = vel_Psmoothed;
      Rvel_[curlevel] = vel_Rsmoothed;
    }
    else                              // PLAIN AGGREGATION
    {
      Pvel_[curlevel] = vel_Ptent;
      Rvel_[curlevel] = vel_Rtent;
    }

    // smooth pressure part
    string preProlongSmoother = preparams->sublist("AMGBS Parameters").get("amgbs: prolongator smoother (pre)","PA-AMG");
    if(preProlongSmoother == "SA-AMG")  // SMOOTHED AGGREGATION
    {
      SA_AMG(A22_[curlevel],pre_Ptent,pre_Rtent,pre_Psmoothed,pre_Rsmoothed);

      Ppre_[curlevel] = pre_Psmoothed;
      Rpre_[curlevel] = pre_Rsmoothed;
    }
    else if(preProlongSmoother == "PG-AMG") // PETROV-GALERKIN
    {
      PG_AMG(A22_[curlevel],pre_Ptent,pre_Rtent,pre_Psmoothed,pre_Rsmoothed);

      Ppre_[curlevel] = pre_Psmoothed;
      Rpre_[curlevel] = pre_Rsmoothed;
    }
    else                                // PLAIN AGGREGATION
    {
      Ppre_[curlevel] = pre_Ptent;
      Rpre_[curlevel] = pre_Rtent;
    }


    if(nVerbose_ > 4) // be verbose
    {
      cout << "Pvel[" << curlevel << "]: " << Pvel_[curlevel]->EpetraMatrix()->NumGlobalRows() << " x " << Pvel_[curlevel]->EpetraMatrix()->NumGlobalCols() << " (" << Pvel_[curlevel]->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
      cout << "Ppre[" << curlevel << "]: " << Ppre_[curlevel]->EpetraMatrix()->NumGlobalRows() << " x " << Ppre_[curlevel]->EpetraMatrix()->NumGlobalCols() << " (" << Ppre_[curlevel]->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;

      cout << "Rvel[" << curlevel << "]: " << Rvel_[curlevel]->EpetraMatrix()->NumGlobalRows() << " x " << Rvel_[curlevel]->EpetraMatrix()->NumGlobalCols() << " (" << Rvel_[curlevel]->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
      cout << "Rpre[" << curlevel << "]: " << Rpre_[curlevel]->EpetraMatrix()->NumGlobalRows() << " x " << Rpre_[curlevel]->EpetraMatrix()->NumGlobalCols() << " (" << Rpre_[curlevel]->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
    }

    /////////////////////////// calc RAP product for next level
    A11_[curlevel+1] = Multiply(*Rvel_[curlevel],*A11_[curlevel],*Pvel_[curlevel]);
    A12_[curlevel+1] = Multiply(*Rvel_[curlevel],*A12_[curlevel],*Ppre_[curlevel]);
    A21_[curlevel+1] = Multiply(*Rpre_[curlevel],*A21_[curlevel],*Pvel_[curlevel]);
    A22_[curlevel+1] = Multiply(*Rpre_[curlevel],*A22_[curlevel],*Ppre_[curlevel]);

    if(nVerbose_ > 4) // be verbose
    {
      cout << "A11[" << curlevel+1 << "]: " << A11_[curlevel+1]->EpetraMatrix()->NumGlobalRows() << " x " << A11_[curlevel+1]->EpetraMatrix()->NumGlobalCols() << " (" << A11_[curlevel+1]->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
      cout << "A12[" << curlevel+1 << "]: " << A12_[curlevel+1]->EpetraMatrix()->NumGlobalRows() << " x " << A12_[curlevel+1]->EpetraMatrix()->NumGlobalCols() << " (" << A12_[curlevel+1]->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
      cout << "A21[" << curlevel+1 << "]: " << A21_[curlevel+1]->EpetraMatrix()->NumGlobalRows() << " x " << A21_[curlevel+1]->EpetraMatrix()->NumGlobalCols() << " (" << A21_[curlevel+1]->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
      cout << "A22[" << curlevel+1 << "]: " << A22_[curlevel+1]->EpetraMatrix()->NumGlobalRows() << " x " << A22_[curlevel+1]->EpetraMatrix()->NumGlobalCols() << " (" << A22_[curlevel+1]->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
    }

    //////////////////// create pre- and postsmoothers
    std::stringstream stream;
    stream << "braess-sarazin: list (level " << curlevel << ")";
    ParameterList& subparams = spparams->sublist("AMGBS Parameters").sublist(stream.str());

    // copy ML Parameters or IFPACK Parameters from FLUID PRESSURE SOLVER block
    if (pressureparams_.isSublist("IFPACK Parameters"))
      subparams.sublist("IFPACK Parameters") = pressureparams_.sublist("IFPACK Parameters");
    else if(pressureparams_.isSublist("ML Parameters"))
      subparams.sublist("ML Parameters") = pressureparams_.sublist("ML Parameters");
    else
      dserror("SaddlePointPreconditioner::Setup: no IFPACK or ML ParameterList found in FLUD PRESSURE SOLVER block -> cannot be!");

    if(curlevel==0)
    {
      subparams.set("pressure correction approx: type",subparams.get("fine: type","IFPACK"));
      if(subparams.isParameter("fine: ifpack type"))  subparams.set("pressure correction approx: ifpack type",subparams.get("fine: ifpack type","ILU"));
      subparams.remove("fine: type",false);
      subparams.remove("fine: ifpack type",false);
    }
    else
    {
      subparams.set("pressure correction approx: type",subparams.get("medium: type","IFPACK"));
      if(subparams.isParameter("medium: ifpack type"))  subparams.set("pressure correction approx: ifpack type",subparams.get("medium: ifpack type","ILU"));
      subparams.remove("medium: type",false);
      subparams.remove("medium: ifpack type",false);
    }

    if(nVerbose_ > 8)
    {
      cout << "Braess-Sarazin smoother (level " << curlevel << ")" << endl << "parameters:" << endl << subparams << endl << endl;
    }

    preS_[curlevel]  = rcp(new BraessSarazin_Smoother(A11_[curlevel],A12_[curlevel],A21_[curlevel],A22_[curlevel],subparams));
    postS_[curlevel] = preS_[curlevel];//rcp(new BraessSarazin_Smoother(A11_[curlevel],A12_[curlevel],A21_[curlevel],A22_[curlevel],subparams));

    //////////////////// prepare variables for next aggregation level
    curvelNS = nextvelNS;
    curpreNS = nextpreNS;

    nlevels_ = curlevel + 1;

    //////////////////// check if aggregation is complete
    if ((A11_[curlevel+1]->EpetraMatrix()->NumGlobalRows() + A22_[curlevel+1]->EpetraMatrix()->NumGlobalRows()) < nmaxcoarsedim)
    {
      if(nVerbose_ > 4) cout << "dim A[" << curlevel+1 << "] < " << nmaxcoarsedim << ". -> end aggregation process" << endl;
      break;
    }
  }

  //////////////////// setup coarsest smoother
  std::stringstream stream;
  stream << "braess-sarazin: list (level " << nlevels_ << ")";
  ParameterList& subparams = spparams->sublist("AMGBS Parameters").sublist(stream.str());

  // copy ML Parameters or IFPACK Parameters from FLUID PRESSURE SOLVER block
  if (pressureparams_.isSublist("IFPACK Parameters"))
    subparams.sublist("IFPACK Parameters") = pressureparams_.sublist("IFPACK Parameters");
  else if(pressureparams_.isSublist("ML Parameters"))
    subparams.sublist("ML Parameters") = pressureparams_.sublist("ML Parameters");
  else
    dserror("SaddlePointPreconditioner::Setup: no IFPACK or ML ParameterList found in FLUD PRESSURE SOLVER block -> cannot be!");

  subparams.set("pressure correction approx: type",subparams.get("coarse: type","IFPACK"));
  if(subparams.isParameter("coarse: ifpack type"))  subparams.set("pressure correction approx: ifpack type",subparams.get("coarse: ifpack type","ILU"));
  subparams.remove("coarse: type",false);
  subparams.remove("coarse: ifpack type",false);

  if(nVerbose_ > 8)
  {
    cout << "Braess-Sarazin smoother (level " << nlevels_ << ")" << endl << "parameters:" << endl << subparams << endl << endl;
  }

  coarsestSmoother_ = rcp(new BraessSarazin_Smoother(A11_[nlevels_],A12_[nlevels_],A21_[nlevels_],A22_[nlevels_],subparams));

  if(nVerbose_ > 2)
  {
    cout << "setup phase complete:" << endl;
    cout << "nlevels/maxlevels: " << nlevels_+1 << "/" << nmaxlevels_+1 << endl;
  }

  MLAPI::Finalize();


#ifdef WRITEOUTSTATISTICS
  if(outfile_)
  {
    fprintf(outfile_,"saddlepointPrecSetupTime %f\tsaddlepointPrecLevels %i\t",ttt.ElapsedTime(),nlevels_);
  }

#ifdef WRITEOUTSYMMETRY
  RCP<SparseMatrix> tmpmtx = rcp(new SparseMatrix(*Ainput_->Merge(),Copy));
  tmpmtx->Add(*Ainput_->Merge(),true,-1.0,1.0);
  fprintf(outfile_,"NormFrobenius %f\t",tmpmtx->NormFrobenius());
#endif
#endif

}



void LINALG::SaddlePointPreconditioner::GetPtent(const Epetra_Map& rowmap, const Epetra_IntVector aggvec, int naggs, ParameterList& List, const Epetra_MultiVector& ThisNS, RCP<Epetra_CrsMatrix>& Ptent, RCP<Epetra_MultiVector>& NextNS, const int domainoffset)
{
  const int nsdim = List.get<int>("null space: dimension",-1);
  if(nsdim <= 0) dserror("null space dimension not given");
  const int mylength = rowmap.NumMyElements();

  //////////////// build a domain map for Ptent
  // find first aggregate on proc
  int firstagg = -1;
  int offset = -1;
  for (int i=0; i<mylength; ++i)
    if (aggvec[i]>=0)
    {
      offset = firstagg = aggvec[i];    // aggregate with agg id "aggvec[i]" is first aggregate on current proc
      break;
    }
  offset *= nsdim;                      // calculate offset with dim of null space
  if (offset < 0) dserror("could not find any aggreagate on proc");

  // calculate gids on coarse grid
  vector<int> coarsegids(naggs*nsdim);
  for (int i=0; i<naggs; ++i)
    for (int j=0; j<nsdim; ++j)
    {
      coarsegids[i*nsdim+j] = offset + domainoffset;
      ++offset;
    }
  Epetra_Map pdomainmap(-1,naggs*nsdim,&coarsegids[0],0,aggvec.Comm()); // this is the coarse grid (domain) map

  ////////////////// loop over aggregates and build ids for dofs
  map<int, vector<int> > aggdofs;
  map<int, vector<int> >::iterator fool;
  for(int i=0; i<naggs; ++i)
  {
    vector<int> gids(0);
    aggdofs.insert(pair<int,vector<int> >(firstagg+i,gids));  // is meant to contain all dof gids for an aggregate
  }
  for(int i=0; i<mylength; ++i)
  {
    if(aggvec[i] < 0) continue;   // this agg doesn't belong to current proc
    vector<int>& gids = aggdofs[aggvec[i]];
    gids.push_back(aggvec.Map().GID(i));  // add dof gid to gids for current aggregate
  }

  //////////////// coarse level nullspace to be filled
  NextNS = rcp(new Epetra_MultiVector(pdomainmap,nsdim,true));
  Epetra_MultiVector& nextns = * NextNS;

  //////////////// create Ptent
  Ptent = rcp(new Epetra_CrsMatrix(Copy, rowmap, nsdim));   // create new Ptent matrix

  // fill Ptent
  // loop over aggregates and extract the appropriate slices of the null space
  // do QR and assemble Q and R to Ptent and NextNS
  for (fool =aggdofs.begin(); fool!=aggdofs.end(); ++fool)
  {
    // extract aggregate-local junk of null space
    const int aggsize = (int) fool->second.size();
    Epetra_SerialDenseMatrix Bagg(aggsize,nsdim);
    for (int i=0; i< aggsize; ++i)
      for(int j=0; j<nsdim; ++j)
        Bagg(i,j) = (*ThisNS(j))[ThisNS.Map().LID(fool->second[i])];

    // Bagg = Q*R
    int m = Bagg.M();
    int n = Bagg.N();
    int lwork = n*10;
    int info = 0;
    int k = min(m,n);
    if(k!=n) dserror("Aggregate too small, fatal!");

    vector<double> work(lwork);
    vector<double> tau(k);
    Epetra_LAPACK lapack;
    lapack.GEQRF(m,n,Bagg.A(),m,&tau[0],&work[0],lwork,&info);
    if (info) dserror ("Lapack dgeqrf returned nonzero");
    if (work[0]>lwork)
    {
      lwork = (int) work[0];
      work.resize(lwork);
    }

    // get R (stored on Bagg) and assemble it into nextns
    int agg_cgid = fool->first*nsdim;
    if(!nextns.Map().MyGID(agg_cgid+domainoffset)) dserror("Missing coarse column id on this proc");
    for (int i=0; i<n; ++i)
      for (int j=i; j<n; ++j)
        (*nextns(j))[nextns.Map().LID(domainoffset+agg_cgid+i)] = Bagg(i,j);

    // get Q and assemble it into Ptent
    lapack.ORGQR(m,n,k,Bagg.A(),m,&tau[0],&work[0],lwork,&info);
    if (info) dserror("Lapack dorgqr returned nonzero");
    for (int i=0; i<aggsize; ++i)
    {
      const int actgrow = fool->second[i];
      for (int j=0; j<nsdim; ++j)
      {
        int actgcol = fool->first*nsdim+j+domainoffset;
        int errone = Ptent->SumIntoGlobalValues(actgrow,1,&Bagg(i,j),&actgcol);
        if (errone>0)
        {
          int errtwo = Ptent->InsertGlobalValues(actgrow,1,&Bagg(i,j),&actgcol);
          if (errtwo<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned negative nonzero");
        }
        else if (errone) dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned negative nonzero");
      }
    } // for (int i=0; i<aggsize; ++i)
  } // for (fool=aggdofs.begin(); fool!=aggdofs.end(); ++fool)
  int err = Ptent->FillComplete(pdomainmap,rowmap);
  if (err) dserror("Epetra_CrsMatrix::FillComplete returned nonzero");
  err = Ptent->OptimizeStorage();
  if (err) dserror("Epetra_CrsMatrix::OptimizeStorage returned nonzero");
}

void LINALG::SaddlePointPreconditioner::GetPtent(const Epetra_Map& rowmap, const Epetra_IntVector aggvec, int naggs, ParameterList& List, const Epetra_MultiVector& ThisNS, RCP<SparseMatrix>& Ptent, RCP<Epetra_MultiVector>& NextNS, const int domainoffset)
{
  const int nsdim = List.get<int>("null space: dimension",-1);
  if(nsdim <= 0) dserror("null space dimension not given");
  const int mylength = rowmap.NumMyElements();

  //////////////// build a domain map for Ptent
  // find first aggregate on proc
  int firstagg = -1;
  int offset = -1;
  for (int i=0; i<mylength; ++i)
    if (aggvec[i]>=0)
    {
      offset = firstagg = aggvec[i];    // aggregate with agg id "aggvec[i]" is first aggregate on current proc
      break;
    }
  offset *= nsdim;                      // calculate offset with dim of null space
  if (offset < 0) dserror("could not find any aggreagate on proc");

  // calculate gids on coarse grid
  vector<int> coarsegids(naggs*nsdim);
  for (int i=0; i<naggs; ++i)
    for (int j=0; j<nsdim; ++j)
    {
      coarsegids[i*nsdim+j] = offset + domainoffset;
      ++offset;
    }
  Epetra_Map pdomainmap(-1,naggs*nsdim,&coarsegids[0],0,aggvec.Comm()); // this is the coarse grid (domain) map

  ////////////////// loop over aggregates and build ids for dofs
  map<int, vector<int> > aggdofs;
  map<int, vector<int> >::iterator fool;
  for(int i=0; i<naggs; ++i)
  {
    vector<int> gids(0);
    aggdofs.insert(pair<int,vector<int> >(firstagg+i,gids));  // is meant to contain all dof gids for an aggregate
  }
  for(int i=0; i<mylength; ++i)
  {
    if(aggvec[i] < 0) continue;   // this agg doesn't belong to current proc
    vector<int>& gids = aggdofs[aggvec[i]];
    gids.push_back(aggvec.Map().GID(i));  // add dof gid to gids for current aggregate
  }

  //////////////// coarse level nullspace to be filled
  NextNS = rcp(new Epetra_MultiVector(pdomainmap,nsdim,true));
  Epetra_MultiVector& nextns = * NextNS;

  //////////////// create Ptent
  Ptent = rcp(new SparseMatrix(rowmap, nsdim));   // create new Ptent matrix

  // fill Ptent
  // loop over aggregates and extract the appropriate slices of the null space
  // do QR and assemble Q and R to Ptent and NextNS
  for (fool =aggdofs.begin(); fool!=aggdofs.end(); ++fool)
  {
    // extract aggregate-local junk of null space
    const int aggsize = (int) fool->second.size();
    Epetra_SerialDenseMatrix Bagg(aggsize,nsdim);
    for (int i=0; i< aggsize; ++i)
      for(int j=0; j<nsdim; ++j)
        Bagg(i,j) = (*ThisNS(j))[ThisNS.Map().LID(fool->second[i])];

    // Bagg = Q*R
    int m = Bagg.M();
    int n = Bagg.N();
    int lwork = n*10;
    int info = 0;
    int k = min(m,n);
    if(k!=n) dserror("Aggregate too small, fatal!");

    vector<double> work(lwork);
    vector<double> tau(k);
    Epetra_LAPACK lapack;
    lapack.GEQRF(m,n,Bagg.A(),m,&tau[0],&work[0],lwork,&info);
    if (info) dserror ("Lapack dgeqrf returned nonzero");
    if (work[0]>lwork)
    {
      lwork = (int) work[0];
      work.resize(lwork);
    }

    // get R (stored on Bagg) and assemble it into nextns
    int agg_cgid = fool->first*nsdim;
    if(!nextns.Map().MyGID(agg_cgid+domainoffset)) dserror("Missing coarse column id on this proc");
    for (int i=0; i<n; ++i)
      for (int j=i; j<n; ++j)
        (*nextns(j))[nextns.Map().LID(domainoffset+agg_cgid+i)] = Bagg(i,j);

    // get Q and assemble it into Ptent
    lapack.ORGQR(m,n,k,Bagg.A(),m,&tau[0],&work[0],lwork,&info);
    if (info) dserror("Lapack dorgqr returned nonzero");
    for (int i=0; i<aggsize; ++i)
    {
      const int actgrow = fool->second[i];
      for (int j=0; j<nsdim; ++j)
      {
        int actgcol = fool->first*nsdim+j+domainoffset;
        int errone = Ptent->EpetraMatrix()->SumIntoGlobalValues(actgrow,1,&Bagg(i,j),&actgcol);
        if (errone>0)
        {
          int errtwo = Ptent->EpetraMatrix()->InsertGlobalValues(actgrow,1,&Bagg(i,j),&actgcol);
          if (errtwo<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned negative nonzero");
        }
        else if (errone) dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned negative nonzero");
      }
    } // for (int i=0; i<aggsize; ++i)
  } // for (fool=aggdofs.begin(); fool!=aggdofs.end(); ++fool)
  int err = Ptent->EpetraMatrix()->FillComplete(pdomainmap,rowmap);
  if (err) dserror("Epetra_CrsMatrix::FillComplete returned nonzero");
  err = Ptent->EpetraMatrix()->OptimizeStorage();
  if (err) dserror("Epetra_CrsMatrix::OptimizeStorage returned nonzero");
}

///////////////////////////////////////////////////////////
int LINALG::SaddlePointPreconditioner::GetGlobalAggregates(SparseMatrix& A, ParameterList& List, const Epetra_MultiVector& ThisNS, Epetra_IntVector& aggrinfo, int& naggregates_local)
{
  int naggregates = GetAggregates(A,List,ThisNS,aggrinfo);

  const Epetra_Comm& comm = A.Comm();
  vector<int> local(comm.NumProc());
  vector<int> global(comm.NumProc());
  for (int i=0; i<comm.NumProc(); ++i) local[i] = 0;  // zero out local vector
  local[comm.MyPID()] = naggregates;                  // fill in local aggregates
  comm.SumAll(&local[0],&global[0],comm.NumProc());   // now all aggregates are in global
  int offset = 0;
  for (int i=0; i<comm.MyPID(); ++i) offset += global[i];
  for (int i=0; i<aggrinfo.MyLength(); ++i)
    if (aggrinfo[i] < naggregates) aggrinfo[i] += offset; // shift "local" agg id to "global" agg id
    else                           aggrinfo[i] = -1;      // set agg info of all non local dofs to -1

  int naggregatesglobal = 0;
  for (int i=0; i<comm.NumProc(); ++i)    // sum up all number of aggregates over all processors
  {
    naggregatesglobal += global[i];
  }

  naggregates_local = naggregates;  // return local number of aggregates for current processor as reference
  return naggregatesglobal;
}

int LINALG::SaddlePointPreconditioner::GetAggregates(SparseMatrix& A, ParameterList& List, const Epetra_MultiVector& ThisNS, Epetra_IntVector& aggrinfo)
{
  if(!A.RowMap().SameAs(aggrinfo.Map())) dserror ("map of aggrinfo must match row map of operator");

  string CoarsenType    = List.get("aggregation: type","Uncoupled");
  double Threshold    = List.get("aggregation: threshold", 0.0);
  int NumPDEEquations   = List.get("PDE equations",1);
  int nsdim         = List.get("null space: dimension", -1);
  if (nsdim==-1)  cout << "dimension of null space not set" << endl;
  int size = A.RowMap().NumMyElements();

  // create ML objects
  ML_Aggregate* agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_KeepInfo(agg_object,1);
  ML_Aggregate_Set_MaxLevels(agg_object,2);
  ML_Aggregate_Set_StartLevel(agg_object,0);
  ML_Aggregate_Set_Threshold(agg_object,Threshold);

  ML_Set_PrintLevel(List.get("ML output", 0));

  // create ML operator
  ML_Operator* ML_Ptent = 0;
  ML_Ptent = ML_Operator_Create(MLAPI::GetML_Comm());

  //if(!thisns) cout << "error: null space is NULL" << endl;
  if (ThisNS.NumVectors() == 0) dserror("error: zero-dimension null space");

  int ns_size = ThisNS.MyLength();

  double* null_vect = 0;
  ML_memory_alloc((void **)(&null_vect), sizeof(double) * ns_size * ThisNS.NumVectors(), "ns");

  int incr = 1;
  for (int v = 0 ; v < ThisNS.NumVectors() ; ++v)
    DCOPY_F77(&ns_size, (double*)ThisNS[v], &incr,
        null_vect + v * ThisNS.MyLength(), &incr);

  ML_Aggregate_Set_NullSpace(agg_object,NumPDEEquations,nsdim,null_vect,size);

  // set coarsening type
  if(CoarsenType == "Uncoupled")
    agg_object->coarsen_scheme = ML_AGGR_UNCOUPLED;
  else if (CoarsenType == "Uncoupled-MIS")
    agg_object->coarsen_scheme = ML_AGGR_HYBRIDUM;
  else if(CoarsenType == "MIS")
  { // needed for MIS, otherwise it sets the number of equations to the null space dimension
    //agg_object->max_levels = -7; // i don't understand this
    agg_object->coarsen_scheme = ML_AGGR_MIS;
  }
  else if(CoarsenType == "METIS")
    agg_object->coarsen_scheme = ML_AGGR_METIS;
  else
  {
    dserror(std::string("error: requested aggregation scheme (" + CoarsenType + ") not recognized"));
  }

  // create ML_Operator for A
  ML_Operator* ML_A = ML_Operator_Create(MLAPI::GetML_Comm());
  ML_Operator_WrapEpetraMatrix(A.EpetraMatrix().get(),ML_A);

  // run coarsening process
  int NextSize = ML_Aggregate_Coarsen(agg_object, ML_A, &ML_Ptent, MLAPI::GetML_Comm());

  int* aggrmap = NULL;
  ML_Aggregate_Get_AggrMap(agg_object,0,&aggrmap);
  if (!aggrmap) dserror("aggr_info not available");

#if 0 // debugging
  fflush(stdout);
  for (int proc=0; proc<A.GetRowMatrix()->Comm().NumProc(); ++proc)
  {
    if (A.GetRowMatrix()->Comm().MyPID()==proc)
    {
      cout << "Proc " << proc << ":" << endl;
      cout << "aggrcount " << aggrcount << endl;
      cout << "NextSize " << NextSize << endl;
      for (int i=0; i<size; ++i)
        cout << "aggrmap[" << i << "] = " << aggrmap[i] << endl;
      fflush(stdout);
    }
    A.GetRowMatrix()->Comm().Barrier();
  }
#endif

  assert (NextSize * nsdim != 0);
  for (int i=0; i<size; ++i) aggrinfo[i] = aggrmap[i];

  ML_Aggregate_Destroy(&agg_object);

  ////////////////////////////////
  // -> i think, Michael forgot this
  // since we're only interested in the aggregates we can free the ML_Operators
  // now valgrind isn't complaining any more
  // but there are still two reachable blocks for Uncoupled coarsening scheme (in ml_qr_fix 15 and 20, called by ML_Aggregate_CoarsenUncoupled in line 629, ml_agg_uncoupled.c)
  // i think it is as ML_qr_fix_setNumDeadNod(numDeadNod); is never called???
  ML_Operator_Destroy(&ML_Ptent);
  ML_Operator_Destroy(&ML_A);
  ML_Ptent = NULL;
  ML_A = NULL;
  ML_qr_fix_Destroy();   // <- ok, this is missing in ML_Aggregate_CoarsenUncoupled in line 629, ml_agg_uncoupled.c

  ML_memory_free((void**)(&null_vect));  // temporary vector with null space data
  null_vect = NULL;
  ////////////////////////////////

  return (NextSize/nsdim);
}

// adapted from MLAPI_Eig.cpp
// if DiagonalScaling == false -> calc max EW of A
// if DiagonalScaling == true  -> calc max EW of D^{-1}A where D is the diagonal of A
// MLAPI has to be initialized before (call MLAPI::Init)
double LINALG::SaddlePointPreconditioner::MaxEigCG(const SparseMatrix& A, const bool DiagonalScaling)
{
  ML_Krylov* kdata = NULL;
  ML_Operator* ML_A = NULL;
  double MaxEigenvalue = 0.0;

  try
  {
    //TEUCHOS_FUNC_TIME_MONITOR("SaddlePointPreconditioner::MaxEigCG");

    // create ML_Operator from SparseMatrix A
    ML_A = ML_Operator_Create(MLAPI::GetML_Comm());
    ML_Operator_WrapEpetraMatrix(A.EpetraMatrix().get(),ML_A);

    kdata = ML_Krylov_Create(MLAPI::GetML_Comm());

    if(DiagonalScaling == false)
      kdata->ML_dont_scale_by_diag = ML_TRUE;
    else
      kdata->ML_dont_scale_by_diag = ML_FALSE;
    ML_Krylov_Set_PrintFreq(kdata,0);
    ML_Krylov_Set_ComputeEigenvalues(kdata);
    ML_Krylov_Set_Amatrix(kdata, ML_A);
    ML_Krylov_Solve(kdata, ML_A->outvec_leng, NULL, NULL);
    MaxEigenvalue = ML_Krylov_Get_MaxEigenvalue(kdata);

    if(MaxEigenvalue == 0.0)  throw std::string("error in MaxEigCG");

    ML_Krylov_Destroy(&kdata);
    ML_Operator_Destroy(&ML_A);
    ML_A = NULL;
    kdata = NULL;

    return MaxEigenvalue;
  }
  catch(std::string str)
  {
    cout << "try to free memory" << endl;
    if(kdata!=NULL) ML_Krylov_Destroy(&kdata); kdata = NULL;
    if(ML_A!=NULL) ML_Operator_Destroy(&ML_A); ML_A = NULL;

    dserror(str);
  }
  return -1.0;
}

// smoothed aggregation (SA-AMG)
void LINALG::SaddlePointPreconditioner::SA_AMG(const RCP<SparseMatrix>& A, const RCP<SparseMatrix>& P_tent, const RCP<SparseMatrix>& R_tent, RCP<SparseMatrix>& P_smoothed, RCP<SparseMatrix>& R_smoothed)
{
  TEUCHOS_FUNC_TIME_MONITOR("SaddlePoint_Preconditioner::SA_AMG");

  ////////////////////
  /*for(int row=0; row < P_tent->EpetraMatrix()->NumMyRows(); row++)
  {
      int numEntries = P_tent->EpetraMatrix()->NumMyEntries(row);
      int indices[numEntries];
      double vals[numEntries];
      int realNumEntries;
      P_tent->EpetraMatrix()->ExtractMyRowCopy(row,numEntries,realNumEntries,&vals[0],&indices[0]);
      cout << "ROW " << row << " numEntries " << numEntries << " rNumEntries " << realNumEntries << " VAL0: " << vals[0] << " INDEX0 " << indices[0] << endl;
      if (numEntries == 0 || realNumEntries == 0)
        cout << "ROW " << row << " has no entry" << endl;
      bool bLineOK = false;
      for(int k=0; k<realNumEntries; k++)
      {
        if(abs(vals[k])>0.0000001)
        {
          bLineOK = true;
          break;
        }
      }
      if(bLineOK == false)
      {
        cout << "ROW " << row << " seems to be a zero row " << endl;
      }
  }*/

  ////////////////////

  double dampingFactor = 1.3333333; // TODO parameter....

  ////////////////// calculate max eigenvalue of diagFinvA (this is a MLAPI call)
  double maxeig = MaxEigCG(*A,true);

  ////////////////// extract diagonal of A
  RCP<Epetra_Vector> diagA = rcp(new Epetra_Vector(A->RowMap(),true));
  A->ExtractDiagonalCopy(*diagA);

  int err = diagA->Reciprocal(*diagA);
  if(err) dserror("SaddlePointPreconditioner::SA_AMG: diagonal entries of A are 0");

#if 0
  // create diagonal matrix
  RCP<SparseMatrix> diagFinv = Teuchos::rcp(new SparseMatrix(A->RowMap(),1,true));
  for (int i=0; i<diagA->Map().NumMyElements(); ++i)
  {
    int gid = diagA->Map().GID(i);
    double val = (*diagA)[i];
    int err = diagFinv->EpetraMatrix()->InsertGlobalValues(gid,1,&val,&gid);
    if (err < 0) dserror(std::string("Epetra_CrsMatrix::InsertGlobalValues returned error code" + err ));
  }
  diagFinv->Complete();

  diagFinv->Scale(dampingFactor/maxeig);  // now we have damping/maxeig(D^{-1} A) * D^{-1}

  diagA = null; // free diagA

  P_smoothed = rcp(new SparseMatrix(P_tent->RowMap(),16,false));   // output = Ptent
  P_smoothed->Add(*P_tent,false,1.0,0.0);

  RCP<SparseMatrix> diagFinvAPtent = Multiply(0,*diagFinv,false,*A,false,*P_tent,false,true); // uses MLMultiply :-)
  P_smoothed->Add(*diagFinvAPtent,false,-1.0,1.0); // P_{SA-AMG} = Ptent - damping/maxeig(D^{-1} A) * D^{-1} A Ptent
  P_smoothed->Complete(P_tent->DomainMap(),P_tent->RangeMap());
#else

  RCP<SparseMatrix> Ascaled = rcp(new SparseMatrix(*A,Copy)); // ok, not the best but just works
  diagA->Scale(dampingFactor/maxeig);
  Ascaled->LeftScale(*diagA);                               // Ascaled = damping/maxeig(D^{-1} A) * D^{-1} * A
  P_smoothed = LINALG::MLMultiply(*Ascaled,*P_tent,false);  // Psmoothed = damping/maxeig(D^{-1} A) * D^{-1} * A * Ptent
  P_smoothed->Add(*P_tent,false,1.0,-1.0);                  // P_smoothed = Ptent - damping/maxeig(D^{-1} A) * D^{-1} * A * Ptent
  P_smoothed->Complete(P_tent->DomainMap(),P_tent->RangeMap());



#endif


  R_smoothed = P_smoothed->Transpose();

}

void LINALG::SaddlePointPreconditioner::PG_AMG(const RCP<SparseMatrix>& A, const RCP<SparseMatrix>& P_tent, const RCP<SparseMatrix>& R_tent, RCP<SparseMatrix>& P_smoothed, RCP<SparseMatrix>& R_smoothed)
{
  TEUCHOS_FUNC_TIME_MONITOR("SaddlePoint_Preconditioner::PG_AMG");


  ////////////////// extract diagonal of A
  RCP<Epetra_Vector> diagA = rcp(new Epetra_Vector(A->RowMap(),true));
  A->ExtractDiagonalCopy(*diagA);
  int err = diagA->Reciprocal(*diagA);
  if(err) dserror("SaddlePointPreconditioner::PG_AMG: diagonal entries of A are 0");


  ///////////////// compute D^{-1}*A
  RCP<SparseMatrix> DinvA = rcp(new SparseMatrix(*A,Copy)); // ok, not the best but just works
  DinvA->LeftScale(*diagA);

  ///////////////// calculate D^{-1}*A*P0
  RCP<SparseMatrix> DinvAP0 = LINALG::MLMultiply(*DinvA,*P_tent,true);

  // TODO: drop values in DinvAP0

  // TODO: compress DinvAP0 -> DinvAP0_subset

  ///////////////// prepare variables for column-based omegas
  int NComputedOmegas = P_tent->EpetraMatrix()->DomainMap().NumGlobalElements();  // number of col-based omegas to be computed (depends on DinvAP0_subset)

  ///////////////// compute D^{-1} * A * D^{-1} * A * P0
  RCP<SparseMatrix> DinvADinvAP0 = LINALG::MLMultiply(*DinvA,*DinvAP0,true);

  // TODO: drop small values in DinvADinvAP0

  ///////////////// define variables for column-based omegas
  std::vector<double> Numerator(NComputedOmegas, 0.0);
  std::vector<double> Denominator(NComputedOmegas, 0.0);
  std::vector<double> ColBasedOmegas(NComputedOmegas, -666.0);

  ///////////////// minimize with respect to the (D^{-1} A)' D^{-1} A norm.
  //
  //               diag( P0' (A' D^{-1}' D^{-1} A) D^{-1} A P0 )
  //  omega = ---------------------------------------------------------
  //           diag( P0' A' D^{-1} (A' D^{-1}' D^{-1} A) D^{-1} A P0 )

#if 1

  // runs only in serial!
  MultiplyAll(DinvAP0,DinvADinvAP0,Numerator);  // -> DinvAP0_subset
  MultiplySelfAll(DinvADinvAP0,Denominator);

  int count_local_zero_denoms = 0;
  for(int i=0; i<Denominator.size(); i++)
    if(Denominator[i]==0.0)
      count_local_zero_denoms++;
  int count_global_zero_denoms = 0;
  count_global_zero_denoms = count_local_zero_denoms;
  //Comm().SumAll(&count_local_zero_denoms,&count_global_zero_denoms,1);
  cout << "zeros: " << count_global_zero_denoms << " out of " << Denominator.size() << endl;


  ///////////////// compute col-based omegas
  int zero_local = 0;   // number of local zeros
  double min_local = DBL_MAX;
  double max_local = DBL_MIN;

  for(int i=0; i<NComputedOmegas; i++)
  {
#ifdef DEBUG
    //if(Denominator[i]==0) cout << "WARNING: Denominator[" << i << "] is zero!!!\n\n" << endl;
#endif
    ColBasedOmegas[i] = Numerator[i]/Denominator[i];  // calculate omegas
    double& val = ColBasedOmegas[i];
    if(val < 0.0)
    {
      val = 0.0;    // zero out negative omegas;
      ++zero_local;
    }
    if(nVerbose_ > 6)
    {
      if(val < min_local) min_local = val;
      if(val > max_local) max_local = val;
    }
  }

  // be verbose
  if(nVerbose_ > 6)
  {
    double min_global = DBL_MAX;
    double max_global = DBL_MIN;
    int zeros_global = 0;
    Comm().MinAll(&min_local,&min_global,1);
    Comm().MaxAll(&max_local,&max_global,1);
    Comm().SumAll(&zero_local,&zeros_global,1);
    cout << "-------------------------------------" << endl;
    cout << "PG-AMG: damping parameters: min="<< min_global << " mx=" << max_global << endl;
    cout << zeros_global << " out of " << NComputedOmegas << endl;
    cout << "-------------------------------------" << endl;
  }
#else

  // new routine for calculating col omega?

  // first: calculate Denominator
  RCP<Epetra_Vector> global_NumeratorVec = rcp(new Epetra_Vector(DinvADinvAP0->DomainMap(),true));
  RCP<Epetra_Vector> global_DenominatorVec = rcp(new Epetra_Vector(DinvADinvAP0->DomainMap(),true));

  // build denominator
  for(int n=0; n<DinvADinvAP0->EpetraMatrix()->NumMyRows(); n++)
  {
    int nnz = DinvADinvAP0->EpetraMatrix()->NumMyEntries(n);
    int indices[nnz];
    double vals[nnz];
    int numEntries;
    DinvADinvAP0->EpetraMatrix()->ExtractMyRowCopy(n,nnz,numEntries,&vals[0],&indices[0]);

    vector<double> resvals(nnz,0.0);
    for(int i=0; i<nnz; i++)
    {
      resvals[i] = (vals[i]*vals[i]);
    }
    global_DenominatorVec->SumIntoMyValues(nnz,&resvals[0],&indices[0]);
  }

  for(int n=0; n<global_DenominatorVec->MyLength(); n++)
    if((*global_DenominatorVec)[n]==0.0)
      cout << "entry is zero" << endl;

  // build numerator
  for(int n=0; n<DinvAP0->EpetraMatrix()->NumMyRows(); n++)
  {
    int nnzl = DinvAP0->EpetraMatrix()->NumMyEntries(n);
    int indicesl[nnzl];
    double valsl[nnzl];
    int numEntries;
    DinvAP0->EpetraMatrix()->ExtractMyRowCopy(n,nnzl,numEntries,&valsl[0],&indicesl[0]);

    int nnzr = DinvADinvAP0->EpetraMatrix()->NumMyEntries(n);
    int indicesr[nnzr];
    double valsr[nnzr];
    DinvADinvAP0->EpetraMatrix()->ExtractMyRowCopy(n,nnzr,numEntries,&valsr[0],&indicesr[0]);

    vector<double> resvals(nnzl,0.0);
    for(int i=0; i<nnzl; i++)
    {
      for(int j=0;j<nnzr; j++)
      {
        if(indicesl[i]==indicesr[j])
        {
          resvals[i] += (valsl[i]*valsr[j]);
        }
      }
    }
    global_NumeratorVec->SumIntoMyValues(nnzl,&resvals[0],&indicesl[0]);
  }

  for(int n=0; n<global_DenominatorVec->MyLength(); n++)
    if((*global_DenominatorVec)[n]==0.0)
      cout << "numerator entry is zero?" << endl;

  if(global_DenominatorVec->MyLength() != global_NumeratorVec->MyLength())
      dserror("vector lengths doesn't match");

  RCP<Epetra_Vector> global_ColBasedOmegaVec = rcp(new Epetra_Vector(global_DenominatorVec->Map(),true));

  // build column based omega
  for(int i=0; i<global_DenominatorVec->MyLength(); i++)
  {
#ifdef DEBUG
    if((*global_DenominatorVec)[i]==0.0) cout << "WARNING: Denominator[" << i << "] is zero!!!\n\n" << endl;
#endif
    (*global_ColBasedOmegaVec)[i] = (*global_NumeratorVec)[i]/(*global_DenominatorVec)[i];
    //if((*global_ColBasedOmegaVec)[i] < 0.0) (*global_ColBasedOmegaVec)[i] = 0.0;
  }

  // build std::vector<double> ColBasedOmegas(NComputedOmegas, -666.0); list
  if(global_ColBasedOmegaVec->ExtractCopy(&ColBasedOmegas[0])!=0) dserror("ExtractCopy failed");
  for(int i=0; i<ColBasedOmegas.size(); i++)
    if(ColBasedOmegas[i]==0.0)
      dserror("ColBasedOmegas is zero");

#endif



  //////////////////// convert column based omegas to row-based omegas
  std::vector<double> RowOmega_local(DinvAP0->RowMap().NumMyElements(), -666.0);  // contains local row omegas of current proc (if not defined -666.0)
  std::vector<int> RowOmegaGrids_local(DinvAP0->RowMap().NumMyElements(), -1);    // contains grids for row omegas of current proc

  ////////////////// calculate max eigenvalue of diagFinvA (this is a MLAPI call)
  double maxeig = MaxEigCG(*A,true); // this is for dirichlet bc's and only for test

  // number of row omegas
  int NRowOmegaSet = 0;

  // loop over processors matrix rows
  for(int row=0; row<DinvAP0->EpetraMatrix()->NumMyRows(); row++)
  {
    // extract local row information for DinvAP0
    int nnz = DinvAP0->EpetraMatrix()->NumMyEntries(row);
    std::vector<int> indices(nnz);   //int indices[nnz];
    std::vector<double> vals(nnz);   //double vals[nnz];
    int numEntries;
    DinvAP0->EpetraMatrix()->ExtractMyRowCopy(row,nnz,numEntries,&vals[0],&indices[0]);

    bool bAtLeastOneDefined = false;

    if(DinvAP0->EpetraMatrix()->GRID(row) != DinvAP0->EpetraMatrix()->IndexBase()-1)
      RowOmegaGrids_local[row] = DinvAP0->EpetraMatrix()->GRID(row); // store grid to local row omega
    else
      cout << "WARNING: there's a row that belongs not to the current proc?? " << endl;

    // execption: DinvAP0 has zero rows
    // should not be, but sometimes may happen if P_tent has zero rows (due to lost nodes in funny MIS aggregation?)
    if(nnz == 0)
    {
      // check if the reason is an empty row in P_tent
      if(P_tent->EpetraMatrix()->NumMyEntries(row) == 0)
      {
        // ok, that is the bad -1 aggregate thing
        //NRowOmegaSet++;
        //RowOmega_local[row] = 0.0;  // per default no damping for this

        // folgende 3 Zeilen bringen nix
        double dampingFactor = 1.3333333; // these are the -1 aggregate dirichlet bc nodes
        bAtLeastOneDefined = true;
        RowOmega_local[row] = dampingFactor/maxeig;
      }
      else
      {
        // very bad :-(
        if(nVerbose_>6) cout << "WARNING: DinvAP0 has a spurious zero row, even though P_tent has nonzero entries! (row lid=" << row << " grid=" << DinvAP0->EpetraMatrix()->GRID(row) << "); RowOmega=0.0" << endl;
#ifdef DEBUG
        // check A
        if(A->EpetraMatrix()->NumMyEntries(row) == 0)
          cout << "A seems to be the bad guy and has a zero row." << endl;
#endif
        //NRowOmegaSet++;
        //RowOmega_local[row] = 0.0;  // per default no damping for this

        double dampingFactor = 1.3333333; // these are the -1 aggregate dirichlet bc nodes
        bAtLeastOneDefined = true;
        RowOmega_local[row] = dampingFactor/maxeig;
      }
    }
    else if(nnz == 1) // 12.04.2010 this is new stuff: special treatment of dirichlet bc's
    {
      // this seems to be a dirichlet bc node
      // the loop below would set RowOmega_local to 1.0, but then
      // we obtain more zero rows in the smoothed prolongator

      // here we set RowOmega_local to the values of SA-AMG???
      double dampingFactor = 1.3333333; // TODO parameter....

      bAtLeastOneDefined = true;
      RowOmega_local[row] = dampingFactor/maxeig;
    }
    else
    {
      // loop over all local column entries of current local row
      for(int j=0; j<nnz; j++)
      {
        // find corresponding col-based omega
        double omega = ColBasedOmegas[indices[j]];
        if(omega != -666.0)
        {
          // col based omega is set
          bAtLeastOneDefined = true;

          if(RowOmega_local[row] == -666.0) RowOmega_local[row] = omega;
          else if(omega < RowOmega_local[row]) RowOmega_local[row] = omega;

          /*double dampingFactor = 1.3333333; // TODO parameter....

          bAtLeastOneDefined = true;
          RowOmega_local[row] = dampingFactor/maxeig;*/
          //bAtLeastOneDefined = true;
          //RowOmega_local[row] = 23.0;

        }
        else
        {
          cout << "WARNING: col-based omega is not set!!!" << endl;
        }
      }
    }

    // check new row omega
    if(bAtLeastOneDefined == true)
    {
      NRowOmegaSet++;
      if(RowOmega_local[row] < 0.) RowOmega_local[row] = 0.;
      //cout << RowOmega_local[row] << endl;
    }
  }


  //////////////// fill RowOmega Epetra vector with local RowOmega_local information
  RCP<Epetra_Vector> RowOmega = rcp(new Epetra_Vector(DinvAP0->RowMap(),true));
  RowOmega->PutScalar(-666.0);

  err = RowOmega->ReplaceGlobalValues(RowOmega_local.size(),&RowOmega_local[0],&RowOmegaGrids_local[0]);
  if(err!=0) dserror("ReplaceGlobalValues returned error %i, one or more indices are not associated with calling proc?",err);

#ifdef DEBUG
  // check if RowOmega has -666. entries -> must not be!
  for(int i=0; i<RowOmega->Map().NumMyElements(); i++)
  {
    if((*RowOmega)[i] == -666.0)
    {
      cout << "proc: " << Comm().MyPID() << " RowOmega LID=" << i << " GID=" << RowOmega->Map().GID(i) << endl;
      dserror("RowOmega still has -666. entires -> must not be!");
    }
  }

  // check if all RowOmega's are set
  if(NRowOmegaSet != RowOmega->Map().NumMyElements())
    dserror("NRowOmegaSet and length of RowOmega doesn't match");

#ifdef WRITEOUTSTATISTICS
  if(outfile_)
  {
    double minVal = 0.0; double maxVal = 0.0;
    RowOmega->MinValue(&minVal); RowOmega->MaxValue(&maxVal);
    fprintf(outfile_,"minRowOmega_%i %f\tmaxRowOmega_%i %f\t",RowOmega->Map().NumGlobalElements(), minVal, RowOmega->Map().NumGlobalElements(), maxVal);
  }
#endif

#endif

  //////////////// compute new prolongator
  /*P_smoothed = LINALG::MLMultiply(*A,*P_tent,false);
  P_smoothed->LeftScale(*diagA);
  P_smoothed->LeftScale(*RowOmega);
  P_smoothed->Add(*P_tent,false,1.0,-0.99);  // TODO -1.0 -> sehr viele nullzeilen wegen Dirichleträndern?
  P_smoothed->Complete(P_tent->DomainMap(),P_tent->RangeMap());*/

  RCP<SparseMatrix> OmegaDinvA = rcp(new SparseMatrix(*A,Copy));
  RCP<Epetra_Vector> diagScaling = rcp(new Epetra_Vector(diagA->Map(),true));
  diagScaling->Multiply(1.0,*RowOmega,*diagA,0.0);
  OmegaDinvA->LeftScale(*diagScaling);
  P_smoothed = LINALG::MLMultiply(*OmegaDinvA,*P_tent,false);
  P_smoothed->Add(*P_tent,false,1.0,-1.0);
  P_smoothed->Complete(P_tent->DomainMap(),P_tent->RangeMap());

#ifdef DEBUG
  vector<int> smoothedzerosgids;
  int smoothedzeros = 0;
  int tentzeros = 0;
  for(int i=0; i<P_smoothed->EpetraMatrix()->NumMyRows(); i++)
  {
    int nnz = P_smoothed->EpetraMatrix()->NumMyEntries(i);
    int indices[nnz];
    double vals[nnz];
    int numEntries;
    P_smoothed->EpetraMatrix()->ExtractMyRowCopy(i,nnz,numEntries,&vals[0],&indices[0]);

    if (nnz==0) // zeros from -1 aggregate
    {
      smoothedzeros++;
      //smoothedzerosgids.push_back(i);
    }
    else
    {
      bool bNonzero = false;
      for(int j=0;j<nnz;j++)
        if(vals[j]!=0.0) bNonzero=true;

      if(bNonzero==false)
      {
        //smoothedzerosgids.push_back(i); // zeros from dirichlet bc's
        smoothedzeros++;
      }
    }

    int nnz2 = P_tent->EpetraMatrix()->NumMyEntries(i);
    int indices2[nnz2];
    double vals2[nnz2];
    P_tent->EpetraMatrix()->ExtractMyRowCopy(i,nnz2,numEntries,&vals2[0],&indices2[0]);

    if (nnz2==0) tentzeros++;
    else
    {
      bool bNonzero = false;
      for(int j=0;j<nnz2;j++)
        if(vals2[j]!=0.0) bNonzero=true;

      if(bNonzero==false) tentzeros++;
    }
  }

  cout << "Psmoothed: " << smoothedzeros << " Ptent: " << tentzeros << endl;

//      std::ofstream fileout;
//      std::stringstream fileoutstream;
//      fileoutstream << "/home/wiesner/fluid/cubits/dc/aggregates/dc/smoothedzeros.txt";
//      fileout.open(fileoutstream.str().c_str(),ios_base::out);
//      for(int i=0; i<smoothedzerosgids.size(); i++)
//        fileout << i << "\t" << smoothedzerosgids[i] << endl;
//      fileout.flush();
//      fileout.close();

#endif

  //////////////// compute restrictor

  //R_smoothed = R_tent;
  //R_smoothed = P_smoothed->Transpose(); // geht aber wird im unsymmetrischen Teil ausgebremst

  // manche werden gar nicht geglättet :-( -> verwende SA-AMG omegas dafür?

  // R = R - R * A * Dinv * omegarfine
  /*R_smoothed = LINALG::MLMultiply(*R_tent,*A,false);
  R_smoothed->RightScale(*RowOmega);
  R_smoothed->RightScale(*diagA);
  R_smoothed->Add(*R_tent,false,1.0,-0.99);
  R_smoothed->Complete(R_tent->DomainMap(),R_tent->RangeMap());*/

  RCP<SparseMatrix> ADinvOmega = rcp(new SparseMatrix(*A,Copy));
  ADinvOmega->RightScale(*diagScaling);
  R_smoothed = LINALG::MLMultiply(*R_tent,*ADinvOmega,false);
  R_smoothed->Add(*R_tent,false,1.0,-1.0);
  R_smoothed->Complete(R_tent->DomainMap(),R_tent->RangeMap());



#if 0
  ///////////////// compute D^{-1}*A^T
  RCP<SparseMatrix> DinvAT = A->Transpose(); // ok, not the best but just works
  DinvAT->LeftScale(*diagA);

  ///////////////// calculate D^{-1}*A^T*P0
  RCP<SparseMatrix> DinvATP0 = LINALG::MLMultiply(*DinvAT,*P_tent,true);

  // TODO: drop values in DinvAP0

  // TODO: compress DinvAP0 -> DinvAP0_subset

  ///////////////// prepare variables for column-based omegas
  int NComputedOmegas_restrictor = P_tent->EpetraMatrix()->DomainMap().NumGlobalElements();  // number of col-based omegas to be computed (depends on DinvAP0_subset)

  ///////////////// compute D^{-1} * A^T * D^{-1} * A^T * P0
  RCP<SparseMatrix> DinvATDinvATP0 = LINALG::MLMultiply(*DinvAT,*DinvATP0,true);

  // TODO: drop small values in DinvADinvAP0

  ///////////////// define variables for column-based omegas
  std::vector<double> Numerator_restrictor(NComputedOmegas_restrictor, 0.0);
  std::vector<double> Denominator_restrictor(NComputedOmegas_restrictor, 0.0);
  std::vector<double> ColBasedOmegas_restrictor(NComputedOmegas_restrictor, -666.0);

  ///////////////// minimize with respect to the (D^{-1} A)' D^{-1} A norm.
  //
  //               diag( R0 (A D^{-1}' D^{-1} A' D^{-1} A' R0' )
  //  omega = ---------------------------------------------------------
  //           diag( R0 A D^{-1} A D^{-1} D^{-1} A' D^{-1} A' R0' )

  MultiplyAll(DinvATP0,DinvATDinvATP0,Numerator_restrictor);  // -> DinvAP0_subset
  MultiplySelfAll(DinvATDinvATP0,Denominator_restrictor);

  ///////////////// compute col-based omegas
  zero_local = 0;   // number of local zeros
  min_local = DBL_MAX;
  max_local = DBL_MIN;

  for(int i=0; i<NComputedOmegas_restrictor; i++)
  {
#ifdef DEBUG
    if(Denominator_restrictor[i]==0) cout << "WARNING: Denominator[" << i << "] is zero!!!\n\n" << endl;
#endif
    ColBasedOmegas_restrictor[i] = Numerator_restrictor[i]/Denominator_restrictor[i];  // calculate omegas
    double& val = ColBasedOmegas_restrictor[i];
    if(val < 0.0)
    {
      val = 0.0;    // zero out negative omegas;
      ++zero_local;
    }
    if(nVerbose_ > 6)
    {
      if(val < min_local) min_local = val;
      if(val > max_local) max_local = val;
    }
  }

  // be verbose
  if(nVerbose_ > 6)
  {
    double min_global = DBL_MAX;
    double max_global = DBL_MIN;
    int zeros_global = 0;
    Comm().MinAll(&min_local,&min_global,1);
    Comm().MaxAll(&max_local,&max_global,1);
    Comm().SumAll(&zero_local,&zeros_global,1);
    cout << "-------------------------------------" << endl;
    cout << "PG-AMG: damping parameters: min="<< min_global << " mx=" << max_global << endl;
    cout << zeros_global << " out of " << NComputedOmegas_restrictor << endl;
    cout << "-------------------------------------" << endl;
  }

  //////////////////// convert column based omegas to row-based omegas
  std::vector<double> RowOmega_restrictor_local(DinvATP0->RowMap().NumMyElements(), -666.0);  // contains local row omegas of current proc (if not defined -666.0)
  std::vector<int> RowOmegaGrids_restrictor_local(DinvATP0->RowMap().NumMyElements(), -1);    // contains grids for row omegas of current proc

  // number of row omegas
  int NRowOmegaSet_restrictor = 0;

  // loop over processors matrix rows
  for(int row=0; row<DinvATP0->EpetraMatrix()->NumMyRows(); row++)
  {
    // extract local row information for DinvAP0
    int nnz = DinvATP0->EpetraMatrix()->NumMyEntries(row);
    int indices[nnz];
    double vals[nnz];
    int numEntries;
    DinvATP0->EpetraMatrix()->ExtractMyRowCopy(row,nnz,numEntries,&vals[0],&indices[0]);

    bool bAtLeastOneDefined = false;

    if(DinvATP0->EpetraMatrix()->GRID(row) != DinvATP0->EpetraMatrix()->IndexBase()-1)
      RowOmegaGrids_restrictor_local[row] = DinvATP0->EpetraMatrix()->GRID(row); // store grid to local row omega
    else
      cout << "WARNING: there's a row that belongs not to the current proc?? " << endl;

    // loop over all local column entries of current local row
    for(int j=0; j<nnz; j++)
    {
      // find corresponding col-based omega
      double omega = ColBasedOmegas_restrictor[indices[j]];
      if(omega != -666.0)
      {
        // col based omega is set
        bAtLeastOneDefined = true;

        if(RowOmega_restrictor_local[row] == -666.0) RowOmega_restrictor_local[row] = omega;
        else if(omega < RowOmega_restrictor_local[row]) RowOmega_restrictor_local[row] = omega;
      }
      else
      {
        cout << "WARNING: col-based omega is not set!!!" << endl;
      }
    }

    // execption: DinvAP0 has zero rows
    // should not be, but sometimes may happen if P_tent has zero rows (due to lost nodes in funny MIS aggregation?)
    if(nnz == 0)
    {
      // check if the reason is an empty row in P_tent
      if(P_tent->EpetraMatrix()->NumMyEntries(row) == 0)
      {
        NRowOmegaSet_restrictor++;
        RowOmega_restrictor_local[row] = 0.0;  // per default no damping for this
      }
      else
      {
        // very bad :-(
        if(nVerbose_>6) cout << "WARNING: DinvATP0 has a spurious zero row, even though P_tent has nonzero entries! (row lid=" << row << " grid=" << DinvAP0->EpetraMatrix()->GRID(row) << "); RowOmega=0.0" << endl;
#if DEBUG
        // check A
        if(A->EpetraMatrix()->NumMyEntries(row) == 0)
          cout << "A seems to be the bad guy and has a zero row." << endl;
#endif
        NRowOmegaSet_restrictor++;
        RowOmega_restrictor_local[row] = 0.0;  // per default no damping for this
      }
    }

    // check new row omega
    if(bAtLeastOneDefined == true)
    {
      NRowOmegaSet_restrictor++;
      if(RowOmega_restrictor_local[row] < 0.) RowOmega_restrictor_local[row] = 0.;
    }
  }


  //////////////// fill RowOmega Epetra vector with local RowOmega_local information
  RCP<Epetra_Vector> RowOmega_restrictor = rcp(new Epetra_Vector(DinvATP0->RowMap(),true));
  RowOmega_restrictor->PutScalar(-666.0);

  err = RowOmega_restrictor->ReplaceGlobalValues(RowOmega_restrictor_local.size(),&RowOmega_restrictor_local[0],&RowOmegaGrids_restrictor_local[0]);
  if(err!=0) dserror("ReplaceGlobalValues returned error %i, one or more indices are not associated with calling proc?",err);

#ifdef DEBUG
  // check if RowOmega has -666. entries -> must not be!
  for(int i=0; i<RowOmega_restrictor->Map().NumMyElements(); i++)
  {
    if((*RowOmega_restrictor)[i] == -666.0)
    {
      cout << "proc: " << Comm().MyPID() << " RowOmega LID=" << i << " GID=" << RowOmega_restrictor->Map().GID(i) << endl;
      //dserror("RowOmega still has -666. entires -> must not be!");
    }
  }
#endif

  RCP<SparseMatrix> Ascaled = rcp(new SparseMatrix(*A,Copy));
  Ascaled->RightScale(*RowOmega_restrictor);
  Ascaled->RightScale(*diagA);
  R_smoothed = LINALG::MLMultiply(*R_tent,*Ascaled,false);
  R_smoothed->Add(*R_tent,false,1.0,-0.99);
  R_smoothed->Complete(R_tent->DomainMap(),R_tent->RangeMap());
#endif
}



void LINALG::SaddlePointPreconditioner::MultiplySelfAll(const RCP<SparseMatrix>& Op,std::vector<double>& Column2Norm)
{
  //if(Column2Norm.size()>0)  dserror("please make sure, that Column2Norm vector is empty");
  int nCols = Op->DomainMap().NumGlobalElements();
  if((int)Column2Norm.capacity() != nCols)
    Column2Norm.resize(nCols,0.0);

  std::vector<double> Column2Norm_local(nCols,0.0);

  for(int n=0; n<Op->EpetraMatrix()->NumMyRows(); n++)
  {
    int nnz = Op->EpetraMatrix()->NumMyEntries(n);
    std::vector<int> indices(nnz); //int indices[nnz];
    std::vector<double> vals(nnz); //double vals[nnz];
    int numEntries;
    Op->EpetraMatrix()->ExtractMyRowCopy(n,nnz,numEntries,&vals[0],&indices[0]);

    for(int i=0; i<nnz; i++)
    {
      Column2Norm_local[indices[i]] += (vals[i]*vals[i]);
    }
  }

  Op->Comm().SumAll(&Column2Norm_local[0],&Column2Norm[0],nCols);
}



void LINALG::SaddlePointPreconditioner::MultiplyAll(const RCP<SparseMatrix>& left, const RCP<SparseMatrix>& right, std::vector<double>& InnerProd)
{
  //if(InnerProd.size()>0)  dserror("please make sure, that InnerProd vector is empty");
  if(!left->DomainMap().SameAs(right->DomainMap())) dserror("domain map of left and right does not match");
  if(!left->RowMap().SameAs(right->RowMap())) dserror("row map of left and right does not match");

  int nCols = left->EpetraMatrix()->DomainMap().NumGlobalElements();

  if((int)InnerProd.capacity() != nCols)
    InnerProd.resize(nCols,0.0);

  std::vector<double> InnerProd_local(nCols,0.0);

  for(int n=0; n<left->EpetraMatrix()->NumMyRows(); n++)
  {
    int nnzl = left->EpetraMatrix()->NumMyEntries(n);
    std::vector<int> indicesl(nnzl); //int indicesl[nnzl];
    std::vector<double> valsl(nnzl); //double valsl[nnzl];
    int numEntries;
    left->EpetraMatrix()->ExtractMyRowCopy(n,nnzl,numEntries,&valsl[0],&indicesl[0]);

    int nnzr = right->EpetraMatrix()->NumMyEntries(n);
    std::vector<int> indicesr(nnzr); //int indicesr[nnzr];
    std::vector<double> valsr(nnzr); //double valsr[nnzr];
    right->EpetraMatrix()->ExtractMyRowCopy(n,nnzr,numEntries,&valsr[0],&indicesr[0]);

    for(int i=0; i<nnzl; i++)
    {
      for(int j=0;j<nnzr; j++)
      {
        if(indicesl[i]==indicesr[j])
        {
          InnerProd_local[indicesl[i]] += (valsl[i]*valsr[j]);
        }
      }
    }
  }

  left->Comm().SumAll(&InnerProd_local[0],&InnerProd[0],nCols);
}


///////////////////////////////////////////////////////////////////
RCP<LINALG::SparseMatrix> LINALG::SaddlePointPreconditioner::Multiply(const SparseMatrix& A, const SparseMatrix& B, const SparseMatrix& C, bool bComplete)
{
  TEUCHOS_FUNC_TIME_MONITOR("SaddlePoint_Preconditioner::Multiply (with MLMultiply)");

  RCP<SparseMatrix> tmp = LINALG::MLMultiply(B,C,true);
  return LINALG::MLMultiply(A,*tmp,bComplete);
}






#endif // CCADISCRET

