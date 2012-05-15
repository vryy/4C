/*
 * braesssarazin_smoother.cpp
 *
 *  Created on: Feb 22, 2010
 *      Author: wiesner
 */



#include "braesssarazin_smoother.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>

#include "Teuchos_TimeMonitor.hpp"

#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_PointRelaxation.h"

#include "ml_MultiLevelPreconditioner.h"

LINALG::BraessSarazin_Smoother::BraessSarazin_Smoother(RCP<const SparseMatrix> A11, RCP<const SparseMatrix> A12, RCP<const SparseMatrix> A21, RCP<const SparseMatrix> A22, const ParameterList& params)
:
Epetra_Operator(),
F_(A11),
G_(A12),
D_(A21),
Z_(A22),
params_(params),
omega_(1.3),
velmap_(A11->RowMap()),
premap_(A22->RowMap())
{
  diagFinv_ = null;


  Setup();
}

// TODO: enable vector splitting here
int LINALG::BraessSarazin_Smoother::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  dserror("overloaded ApplyInverse not supported. use ApplyInverse with splitted vectors instead.");

  return 0;
}

int LINALG::BraessSarazin_Smoother::ApplyInverse(const Epetra_MultiVector& velrhs, const Epetra_MultiVector& prerhs, Epetra_MultiVector& velsol, Epetra_MultiVector& presol, int level) const
{
  TEUCHOS_FUNC_TIME_MONITOR("BraessSarazin_Smoother::ApplyInverse");

  ////////////////// define some variables
  RCP<Epetra_Vector> velres = rcp(new Epetra_Vector(F_->RowMap(),true));
  RCP<Epetra_Vector> preres = rcp(new Epetra_Vector(Z_->RowMap(),true));

  RCP<Epetra_Vector> vtemp = rcp(new Epetra_Vector(F_->RowMap(),true));
  RCP<Epetra_Vector> ptemp = rcp(new Epetra_Vector(Z_->RowMap(),true));

  for(int k = 0; k<nSweeps_; k++)
  {
    ////////////////// 1) calculate residual
    //if(level==0)
    {
    F_->Apply(velsol,*vtemp);
    G_->Apply(presol,*velres);
    velres->Update(1.0,*vtemp,1.0);  // velres = + F vsol + G presol
    velres->Update(1.0,velrhs,-1.0); // velres = velrhs - F vsol - G presol

    D_->Apply(velsol,*ptemp);
    Z_->Apply(presol,*preres);
    preres->Update(1.0,*ptemp,1.0); // preres = + D vsol + Z presol
    preres->Update(1.0,prerhs,-1.0); // preres = prerhs - D vsol - Z presol
    }
    //else
    /*{
      velres->Update(1.0,velrhs,0.0);
      preres->Update(1.0,prerhs,0.0);
    }*/

    ////////////////// 2) solve for pressure update q = (D Fhatinv G)^{-1} (D Fhatinv velres - omega preres)
    ////////////////// 2.1) calculate rhs
    vtemp->PutScalar(0.0);
    RCP<Epetra_Vector> qrhs = rcp(new Epetra_Vector(Z_->RowMap(),true));
    vtemp->Multiply(1.0,*diagFinv_,*velres,0.0);
    D_->Apply(*vtemp,*qrhs);
    qrhs->Update(omega_,*preres,-1.0);  // qrhs = -D * Fhatinv * velres + omega * preres

    ////////////////// 2.2) "solve" pressure correction equation for pressure update
    RCP<Epetra_Vector> q = rcp(new Epetra_Vector(Z_->RowMap(),true)); // vector for "solution" (pressure update)
    Pp_->ApplyInverse(*qrhs,*q);    // use IFPACK

    ////////////////// 3) update step
    vtemp->PutScalar(0.0);
    G_->Apply(*q,*vtemp);
    vtemp->Update(1.0,*velres,-1.0);    // velres - G * q

    RCP<Epetra_Vector> vx = rcp(new Epetra_Vector(diagFinv_->Map(),true));
    vx->Multiply(1/omega_,*diagFinv_,*vtemp,0.0); // vx = 1/omega_ diagFinv * (velres - G * q)

    ////////////////// adapt solution v: v = v + [vx;q]

    // velocity part: velsol = velsol + vx
    velsol.Update(1.0,*vx,1.0);

    // pressure part: presol = presol + q
    presol.Update(1.0,*q,1.0);
  }

  return 0;
}


bool LINALG::BraessSarazin_Smoother::Setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("BraessSarazin_Smoother::Setup");

  //cout << DRT::Problem::Instance()->FluidPressureSolverParams() << endl;

  nSweeps_ = params_.get("braess-sarazin: sweeps",1);   // number of Braess-Sarazin sweeps (default)
  omega_   = params_.get("braess-sarazin: damping factor",1.3); // damping factor

  ///////////////////////// extract diagonal of F
#if 0
  Epetra_Vector diagvec(F_->RowMap(),true);
  F_->ExtractDiagonalCopy(diagvec);
  int err = diagvec.Reciprocal(diagvec);
  if (err) dserror("diagonal entries = 0");

  // create diagonal matrix
  diagFinv_ = rcp(new SparseMatrix(F_->RowMap(),1));
  for (int i=0; i<diagvec.Map().NumMyElements(); ++i)
  {
    int gid = diagvec.Map().GID(i);
    double val = diagvec[i];
    int err = diagFinv_->EpetraMatrix()->InsertGlobalValues(gid,1,&val,&gid);
    if (err < 0) dserror(std::string("Epetra_CrsMatrix::InsertGlobalValues returned error code" + err ));
  }

  diagFinv_->Complete();

  ////////////////////////// calculate Schur complement
  S_ = Multiply(*rcp_const_cast<SparseMatrix>(D_),false,*rcp_const_cast<SparseMatrix>(diagFinv_),false,*rcp_const_cast<SparseMatrix>(G_),false,false);
  S_->Add(*Z_,false,omega_,-1.0); // S = omega Z - S
  S_->Complete();
#else
  diagFinv_ = rcp(new Epetra_Vector(F_->RowMap(),true));
  F_->ExtractDiagonalCopy(*diagFinv_);
  int err = diagFinv_->Reciprocal(*diagFinv_);
  if (err) dserror("diagonal entries = 0");

  RCP<SparseMatrix> Gscaled = rcp(new SparseMatrix(*G_,Copy)); // ok, not the best but just works
  Gscaled->LeftScale(*diagFinv_);                               // Ascaled = damping/maxeig(D^{-1} A) * D^{-1} * A
  S_ = LINALG::MLMultiply(*D_,*Gscaled,false);
  S_->Add(*Z_,false,omega_,-1.0);
  S_->Complete();

#endif

  /////////////////////// allocate preconditioner/solver for pressure correction equation

  // standard: IFPACK ILU
  string type = params_.get("pressure correction approx: type","ILU");
  if(type=="ILU")
  {
    ParameterList ifpackParams;
    if(params_.isSublist("IFPACK Parameters"))
      ifpackParams = params_.sublist("IFPACK Parameters");
    else
    { // standard parameters for ifpack
      ifpackParams.set("amesos: solver type", "Amesos_Klu");
      ifpackParams.set("fact: drop tolerance",0);
      ifpackParams.set("fact: level-of-fill",0);
      ifpackParams.set("relaxation: damping factor",0);
      ifpackParams.set("schwarz: combine mode","Add");
      ifpackParams.set("schwarz: reordering type","rcm");
    }

    Ifpack factory;
    Ifpack_Preconditioner* prec = factory.Create(type,S_->EpetraMatrix().get(),params_.get("Ifpack overlap",0));
    prec->SetParameters(ifpackParams);
    prec->Initialize();
    prec->Compute();
    Pp_ = rcp(prec);
    /*cout << "chosen IFPACK params" << endl;
    cout << ifpackParams << endl;
    cout << "out of parameters" << endl;
    cout << params_ << endl;*/
  }
  else if(type=="Jacobi" || type=="Gauss-Seidel" || type=="symmetric Gauss-Seidel")
  {
    Ifpack factory;
    Ifpack_Preconditioner* prec = factory.Create("point relaxation",S_->EpetraMatrix().get(),params_.get("Ifpack overlap",0));
    prec->SetParameters(params_.sublist("IFPACK Parameters"));
    prec->Initialize();
    prec->Compute();
    Pp_ = rcp(prec);
    //cout << *prec << endl;
  }
  else if(type=="Jacobi stand-alone" || type=="Gauss-Seidel stand-alone" || type=="symmetric Gauss-Seidel stand-alone")
  {
    Ifpack factory;
    Ifpack_Preconditioner* prec = factory.Create("point relaxation stand-alone",S_->EpetraMatrix().get(),params_.get("Ifpack overlap",0));
    prec->SetParameters(params_.sublist("IFPACK Parameters"));
    prec->Initialize();
    prec->Compute();
    Pp_ = rcp(prec);
    //cout << *prec << endl;
  }
  else if(type=="KLU")
  {
    ParameterList amesosParams;
    amesosParams.set("amesos: solver type","Amesos_Klu");
    RCP<Ifpack_Amesos> prec = rcp(new Ifpack_Amesos(S_->EpetraMatrix().get()));
    prec->SetParameters(amesosParams);
    prec->Initialize();
    prec->Compute();
    Pp_ = prec;
  }
  else if(type=="Umfpack")
  {
    ParameterList amesosParams;
    amesosParams.set("amesos: solver type","Amesos_Umfpack");
    RCP<Ifpack_Amesos> prec = rcp(new Ifpack_Amesos(S_->EpetraMatrix().get()));
    prec->SetParameters(amesosParams);
    prec->Initialize();
    prec->Compute();
    Pp_ = prec;
  }
  else if(type=="ML")
  {
    /*if(params_.isSublist("ML Parameters"))
    {
      Pp_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*S_->EpetraMatrix(),params_.sublist("ML Paramters")));
    }
    else
      dserror("ML parameters for pressure correction missing!);*/
    dserror("ML is not supported for smoothing any more");
  }
  else
    dserror("we need at least a ML Parameters or an IFPACK Parameters sublist for the approximative solution of the pressure correction equation");

  return true;
}

RCP<LINALG::SparseMatrix> LINALG::BraessSarazin_Smoother::Multiply(const SparseMatrix& A, bool transA, const SparseMatrix& B, bool transB, const SparseMatrix& C, bool transC, bool bComplete)
{
  // it's ok to use MLMultiply here (only quadratic matrices with the same map)
#if 0
  RCP<SparseMatrix> tmp = LINALG::Multiply(B,transB,C,transC,true);
  return LINALG::Multiply(A,transA,*tmp,false,bComplete);
#else
  // use Michael's MLMultiply
  if(transA == true || transB == true || transC == true)
    dserror("up to now we don't support matri-matrix multiplication with transposed flag on");

  RCP<SparseMatrix> tmp = LINALG::MLMultiply(B,C,true);
  return LINALG::MLMultiply(A,*tmp,bComplete);
#endif
}

