/*
 * transfer_operator_saamg.cpp
 *
 *  Created on: Apr 20, 2010
 *      Author: wiesner
 */


#include "transfer_operator.H"
#include "transfer_operator_tentative.H"
#include "transfer_operator_saamg.H"

// includes for MLAPI functions
#include "ml_common.h"
#include "ml_include.h"
#include "ml_aggregate.h"
#include "ml_agg_METIS.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "MLAPI.h"


LINALG::SAAMGTransferOperator::SAAMGTransferOperator(const RCP<SparseMatrix>& A, FILE* outfile) :
  TentativeTransferOperator(A, outfile)
{

}

RCP<Epetra_MultiVector>  LINALG::SAAMGTransferOperator::buildTransferOperators(const RCP<Epetra_IntVector> aggs, int naggs_local, Teuchos::ParameterList& params, const RCP<Epetra_MultiVector>& ThisNS, const int domainoffset)
{
  TEUCHOS_FUNC_TIME_MONITOR("SAAMGTransferOperator::buildTransferOperators");

  ////////////// define dummy variable for next nullspace
  RCP<Epetra_MultiVector> NextNS = null;
  RCP<SparseMatrix> prolongator_tent = null;
  RCP<SparseMatrix> restrictor_tent = null;

  ////////////// build tentative prolongator
  GetPtent(A_->RowMap(),*aggs,naggs_local,params,*ThisNS,prolongator_tent,NextNS,domainoffset);

  ////////////// build tentative restrictor
  restrictor_tent = prolongator_tent->Transpose();

  ////////////// constant part for SA-AMG damping parameter
  double dampingFactor = 1.3333333;

  ////////////////// calculate max eigenvalue of diagFinvA (this is a MLAPI call)
  double maxeig = MaxEigCG(*A_,true);

  ////////////////// extract diagonal of A
  RCP<Epetra_Vector> diagA = rcp(new Epetra_Vector(A_->RowMap(),true));
  A_->ExtractDiagonalCopy(*diagA);

  int err = diagA->Reciprocal(*diagA);
  if(err) dserror("SaddlePointPreconditioner::SA_AMG: diagonal entries of A are 0");

  /////////////////// setup smoothed aggregation prolongator
  RCP<SparseMatrix> Ascaled = rcp(new SparseMatrix(*A_,Copy)); // ok, not the best but just works
  diagA->Scale(dampingFactor/maxeig);
  Ascaled->LeftScale(*diagA);                               // Ascaled = damping/maxeig(D^{-1} A) * D^{-1} * A
  prolongator_ = LINALG::MLMultiply(*Ascaled,*prolongator_tent,false);  // Psmoothed = damping/maxeig(D^{-1} A) * D^{-1} * A * Ptent
  prolongator_->Add(*prolongator_tent,false,1.0,-1.0);                  // P_smoothed = Ptent - damping/maxeig(D^{-1} A) * D^{-1} * A * Ptent
  prolongator_->Complete(prolongator_tent->DomainMap(),prolongator_tent->RangeMap());

  /////////////////// setup restrictor
  restrictor_ = prolongator_->Transpose();

  ////////////// return next level nullspace
  return NextNS;
}

/////////////////////////////////////////////////////////////////////
double LINALG::SAAMGTransferOperator::MaxEigCG(const SparseMatrix& A, const bool DiagonalScaling)
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

