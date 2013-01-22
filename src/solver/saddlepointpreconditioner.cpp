/*
 * saddlepointpreconditioner.cpp
 *
 *  Created on: Feb 16, 2010
 *      Author: wiesner
 */



#define WRITEOUTSTATISTICS
#undef WRITEOUTSYMMETRY // write out || A - A^T ||_F
#define WRITEOUTAGGREGATES

//#define ANALYSIS

#include <Epetra_LocalMap.h>

#include "../linalg/linalg_sparsematrix.H"

#include "saddlepointpreconditioner.H"
#include "amg/transfer_operator.H"
#include "amg/transfer_operator_tentative.H"
#include "amg/transfer_operator_saamg.H"
#include "amg/transfer_operator_pgamg.H"
#include "amg/braesssarazin_smoother.H"
#include "amg/aggregation_method_uncoupled.H"
#include "amg/aggregation_method_ml.H"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TimeMonitor.hpp"


// includes for MLAPI functions
#include "ml_common.h"
#include "ml_include.h"
#include "ml_aggregate.h"
#include "ml_agg_METIS.h"
#include "Teuchos_RCP.hpp"
#include "MLAPI_Error.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI.h"

#include "float.h" // for DBL_MAX and DBL_MIN

LINALG::SaddlePointPreconditioner::SaddlePointPreconditioner(RCP<Epetra_Operator> A, const Teuchos::ParameterList& params, FILE* outfile)
:
Epetra_Operator(),
label_("LINALG::SaddlePointPreconditioner"),
params_(params),
outfile_(outfile),
nVerbose_(0)
{
  Setup(A,params);
}

LINALG::SaddlePointPreconditioner::~SaddlePointPreconditioner()
{

}


int LINALG::SaddlePointPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  TEUCHOS_FUNC_TIME_MONITOR("SaddlePointPreconditioner::ApplyInverse");

  // VCycle

  // note: Aztec might pass X and Y as physically identical objects,
  // so we better deep copy here

  RCP<LINALG::ANA::Vector> Xv = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(0),false));
  RCP<LINALG::ANA::Vector> Xp = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(1),false));

  RCP<LINALG::ANA::Vector> Yv = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(0),/*false*/true));
  RCP<LINALG::ANA::Vector> Yp = Teuchos::rcp(new LINALG::ANA::Vector(*mmex_.Map(1),/*false*/true));

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

#ifdef ANALYSIS
#ifdef DEBUG

  RCP<Epetra_Vector> velresXX = Teuchos::rcp(new Epetra_Vector(Yvel.Map(),true));
  RCP<Epetra_Vector> preresXX = Teuchos::rcp(new Epetra_Vector(Ypre.Map(),true));

  RCP<Epetra_Vector> vtempXX = Teuchos::rcp(new Epetra_Vector(Yvel.Map(),true));
  RCP<Epetra_Vector> ptempXX = Teuchos::rcp(new Epetra_Vector(Ypre.Map(),true));

  A11_[level]->Apply(Yvel,*vtempXX);
  A12_[level]->Apply(Ypre,*velresXX);
  velresXX->Update(1.0,*vtempXX,1.0);  // velres = + F Zvel + G Zpre
  velresXX->Update(1.0,Xvel,-1.0); // velres = Xvel - F Zvel - G Zpre

  A21_[level]->Apply(Yvel,*ptempXX);
  A22_[level]->Apply(Ypre,*preresXX);
  preresXX->Update(1.0,*ptempXX,1.0); // preres = + D Zvel + Z Zpre
  preresXX->Update(1.0,Xpre,-1.0); // preres = Xpre - D Zvel - Z Zpre

  if(A11_[level]->Comm().NumProc() == 1)
  {
    // runs only in serial
    FILE* fp = fopen("analysis.m","w");
    int length = A11_[level]->RowMap().NumGlobalPoints();
    length += A22_[level]->RowMap().NumGlobalPoints();
    fprintf(fp,"A=sparse(%d,%d);\n",length,length);


    // loop over local rows of DinvAP0
    for(int row=0; row<A11_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A11_[level]->EpetraMatrix()->GRID(row);
      int nnz = A11_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A11_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

      numEntries = 0;
      nnz = A12_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      indices.clear(); indices.reserve(nnz);
      vals.clear(); vals.reserve(nnz);
      err = A12_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    for(int row=0; row<A22_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A22_[level]->EpetraMatrix()->GRID(row);
      int nnz = A22_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A22_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

      numEntries = 0;
      nnz = A21_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      indices.clear(); indices.reserve(nnz);
      vals.clear(); vals.reserve(nnz);
      err = A21_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    fprintf(fp,"rhs = zeros(%d,1);\n",length);

    double * MyGlobalElements1 = velresXX->Values();
    for(int i=0; i<velresXX->Map().NumMyElements(); i++)
    {
      int grid = velresXX->Map().GID(i);
      double vv = MyGlobalElements1[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,vv);
    }
    double * MyGlobalElements2 = preresXX->Values();
    for(int i=0; i<preresXX->Map().NumMyElements(); i++)
    {
      int grid = preresXX->Map().GID(i);
      double pp = MyGlobalElements2[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,pp);
    }

    fclose(fp);

    cout << "** BEFORE presmoothing (level " << level << ") **" << endl << "Press y when you're done" << endl;
    char instring[100];
    scanf("%s",instring);
  }
#endif
#endif

  // vectors for presmoothed solution
  RCP<Epetra_MultiVector> Zvel = Teuchos::rcp(new Epetra_MultiVector(Yvel.Map(),1,true));
  RCP<Epetra_MultiVector> Zpre = Teuchos::rcp(new Epetra_MultiVector(Ypre.Map(),1,true));

  Zvel->Update(1.0,Yvel,0.0);
  Zpre->Update(1.0,Ypre,0.0);

  // presmoothing
  if(bPresmoothing_) preS_[level]->ApplyInverse(Xvel,Xpre,*Zvel,*Zpre,level);  // rhs X is fix, initial solution Z = 0 (per definition, see above)
                                                      // note: ApplyInverse expects the "solution" and no solution increment "Delta Z"


  // calculate residual (fine grid)
  RCP<Epetra_Vector> velres = Teuchos::rcp(new Epetra_Vector(Yvel.Map(),true));
  RCP<Epetra_Vector> preres = Teuchos::rcp(new Epetra_Vector(Ypre.Map(),true));

  RCP<Epetra_Vector> vtemp = Teuchos::rcp(new Epetra_Vector(Yvel.Map(),true));
  RCP<Epetra_Vector> ptemp = Teuchos::rcp(new Epetra_Vector(Ypre.Map(),true));

  A11_[level]->Apply(*Zvel,*vtemp);
  A12_[level]->Apply(*Zpre,*velres);
  velres->Update(1.0,*vtemp,1.0);  // velres = + F Zvel + G Zpre
  velres->Update(1.0,Xvel,-1.0); // velres = Xvel - F Zvel - G Zpre

  A21_[level]->Apply(*Zvel,*ptemp);
  A22_[level]->Apply(*Zpre,*preres);
  preres->Update(1.0,*ptemp,1.0); // preres = + D Zvel + Z Zpre
  preres->Update(1.0,Xpre,-1.0); // preres = Xpre - D Zvel - Z Zpre

#ifdef ANALYSIS
#ifdef DEBUG
  if(A11_[level]->Comm().NumProc() == 1)
  {
    // runs only in serial
    FILE* fp = fopen("analysis.m","w");
    int length = A11_[level]->RowMap().NumGlobalPoints();
    length += A22_[level]->RowMap().NumGlobalPoints();
    fprintf(fp,"A=sparse(%d,%d);\n",length,length);


    // loop over local rows of DinvAP0
    for(int row=0; row<A11_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A11_[level]->EpetraMatrix()->GRID(row);
      int nnz = A11_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A11_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

      numEntries = 0;
      nnz = A12_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      indices.clear(); indices.reserve(nnz);
      vals.clear(); vals.reserve(nnz);
      err = A12_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    for(int row=0; row<A22_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A22_[level]->EpetraMatrix()->GRID(row);
      int nnz = A22_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A22_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

      numEntries = 0;
      nnz = A21_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      indices.clear(); indices.reserve(nnz);
      vals.clear(); vals.reserve(nnz);
      err = A21_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    fprintf(fp,"rhs = zeros(%d,1);\n",length);

    double * MyGlobalElements1 = velres->Values();
    for(int i=0; i<velres->Map().NumMyElements(); i++)
    {
      int grid = velres->Map().GID(i);
      double vv = MyGlobalElements1[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,vv);
    }
    double * MyGlobalElements2 = preres->Values();
    for(int i=0; i<preres->Map().NumMyElements(); i++)
    {
      int grid = preres->Map().GID(i);
      double pp = MyGlobalElements2[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,pp);
    }

    fclose(fp);

    cout << "** AFTER presmoothing (level " << level << ") **" << endl << "Press y when you're done" << endl;
    char instring[100];
    scanf("%s",instring);
  }
#endif
#endif

  // calculate coarse residual
  RCP<Epetra_Vector> velres_coarse = Teuchos::rcp(new Epetra_Vector(Tvel_[level]->R().RowMap(),true));
  RCP<Epetra_Vector> preres_coarse = Teuchos::rcp(new Epetra_Vector(Tpre_[level]->R().RowMap(),true));
  Tvel_[level]->R().Apply(*velres,*velres_coarse);
  Tpre_[level]->R().Apply(*preres,*preres_coarse);

  // define vector for coarse level solution
  RCP<Epetra_Vector> velsol_coarse = Teuchos::rcp(new Epetra_Vector(A11_[level+1]->RowMap(),true));
  RCP<Epetra_Vector> presol_coarse = Teuchos::rcp(new Epetra_Vector(A22_[level+1]->RowMap(),true));

  // call Vcycle recursively
  VCycle(*velres_coarse,*preres_coarse,*velsol_coarse,*presol_coarse,level+1);

  // define vectors for prolongated solution
  RCP<Epetra_Vector> velsol_prolongated = Teuchos::rcp(new Epetra_Vector(A11_[level]->RowMap(),true));
  RCP<Epetra_Vector> presol_prolongated = Teuchos::rcp(new Epetra_Vector(A22_[level]->RowMap(),true));

  // prolongate solution
  Tvel_[level]->P().Apply(*velsol_coarse,*velsol_prolongated);
  Tpre_[level]->P().Apply(*presol_coarse,*presol_prolongated);

  // update solution Zvel and Zpre for postsmoother
  Zvel->Update(1.0,*velsol_prolongated,1.0);
  Zpre->Update(1.0,*presol_prolongated,1.0);

#ifdef ANALYSIS
#ifdef DEBUG

  // calculate residual
  velres->PutScalar(0.0);
  preres->PutScalar(0.0);

  vtemp->PutScalar(0.0);
  ptemp->PutScalar(0.0);

  A11_[level]->Apply(*Zvel,*vtemp);
  A12_[level]->Apply(*Zpre,*velres);
  velres->Update(1.0,*vtemp,1.0);  // velres = + F Zvel + G Zpre
  velres->Update(1.0,Xvel,-1.0); // velres = Xvel - F Zvel - G Zpre

  A21_[level]->Apply(*Zvel,*ptemp);
  A22_[level]->Apply(*Zpre,*preres);
  preres->Update(1.0,*ptemp,1.0); // preres = + D Zvel + Z Zpre
  preres->Update(1.0,Xpre,-1.0); // preres = Xpre - D Zvel - Z Zpre

  if(A11_[level]->Comm().NumProc() == 1)
  {
    // runs only in serial
    FILE* fp = fopen("analysis.m","w");
    int length = A11_[level]->RowMap().NumGlobalPoints();
    length += A22_[level]->RowMap().NumGlobalPoints();
    fprintf(fp,"A=sparse(%d,%d);\n",length,length);


    // loop over local rows of DinvAP0
    for(int row=0; row<A11_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A11_[level]->EpetraMatrix()->GRID(row);
      int nnz = A11_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A11_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

      numEntries = 0;
      nnz = A12_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      indices.clear(); indices.reserve(nnz);
      vals.clear(); vals.reserve(nnz);
      err = A12_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    for(int row=0; row<A22_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A22_[level]->EpetraMatrix()->GRID(row);
      int nnz = A22_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A22_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

      numEntries = 0;
      nnz = A21_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      indices.clear(); indices.reserve(nnz);
      vals.clear(); vals.reserve(nnz);
      err = A21_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    fprintf(fp,"rhs = zeros(%d,1);\n",length);

    double * MyGlobalElements1 = velres->Values();
    for(int i=0; i<velres->Map().NumMyElements(); i++)
    {
      int grid = velres->Map().GID(i);
      double vv = MyGlobalElements1[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,vv);
    }
    double * MyGlobalElements2 = preres->Values();
    for(int i=0; i<preres->Map().NumMyElements(); i++)
    {
      int grid = preres->Map().GID(i);
      double pp = MyGlobalElements2[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,pp);
    }

    fclose(fp);

    cout << "** AFTER coarse grid correction (level " << level << ") **" << endl << "Press y when you're done" << endl;
    char instring[100];
    scanf("%s",instring);
  }
#endif
#endif


  // postsmoothing
  if (bPostsmoothing_) postS_[level]->ApplyInverse(Xvel,Xpre,*Zvel,*Zpre,level); // rhs the same as for presmoothing, but better initial solution (Z)


#ifdef ANALYSIS
#ifdef DEBUG

  // calculate residual
  velres->PutScalar(0.0);
  preres->PutScalar(0.0);

  vtemp->PutScalar(0.0);
  ptemp->PutScalar(0.0);

  A11_[level]->Apply(*Zvel,*vtemp);
  A12_[level]->Apply(*Zpre,*velres);
  velres->Update(1.0,*vtemp,1.0);  // velres = + F Zvel + G Zpre
  velres->Update(1.0,Xvel,-1.0); // velres = Xvel - F Zvel - G Zpre

  A21_[level]->Apply(*Zvel,*ptemp);
  A22_[level]->Apply(*Zpre,*preres);
  preres->Update(1.0,*ptemp,1.0); // preres = + D Zvel + Z Zpre
  preres->Update(1.0,Xpre,-1.0); // preres = Xpre - D Zvel - Z Zpre

  if(A11_[level]->Comm().NumProc() == 1)
  {
    // runs only in serial
    FILE* fp = fopen("analysis.m","w");
    int length = A11_[level]->RowMap().NumGlobalPoints();
    length += A22_[level]->RowMap().NumGlobalPoints();
    fprintf(fp,"A=sparse(%d,%d);\n",length,length);


    // loop over local rows of DinvAP0
    for(int row=0; row<A11_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A11_[level]->EpetraMatrix()->GRID(row);
      int nnz = A11_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A11_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

      numEntries = 0;
      nnz = A12_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      indices.clear(); indices.reserve(nnz);
      vals.clear(); vals.reserve(nnz);
      err = A12_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    for(int row=0; row<A22_[level]->EpetraMatrix()->NumMyRows(); row++)
    {
      //////////////// extract global information for local row in DinvAP0
      // we need global column ids
      int grid = A22_[level]->EpetraMatrix()->GRID(row);
      int nnz = A22_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      std::vector<int> indices(nnz);
      std::vector<double> vals(nnz);
      int numEntries;
      int err = A22_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

      numEntries = 0;
      nnz = A21_[level]->EpetraMatrix()->NumGlobalEntries(grid);
      indices.clear(); indices.reserve(nnz);
      vals.clear(); vals.reserve(nnz);
      err = A21_[level]->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
      for(int col=0; col<numEntries; col++)
        fprintf(fp,"A(%d,%d)=%25.16e;\n",grid+1,indices[col]+1,vals[col]);

    }

    fprintf(fp,"rhs = zeros(%d,1);\n",length);

    double * MyGlobalElements1 = velres->Values();
    for(int i=0; i<velres->Map().NumMyElements(); i++)
    {
      int grid = velres->Map().GID(i);
      double vv = MyGlobalElements1[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,vv);
    }
    double * MyGlobalElements2 = preres->Values();
    for(int i=0; i<preres->Map().NumMyElements(); i++)
    {
      int grid = preres->Map().GID(i);
      double pp = MyGlobalElements2[i];
      fprintf(fp,"rhs(%d) = %25.16e;\n",grid+1,pp);
    }

    fclose(fp);

    cout << "** AFTER post smoothing (level " << level << ") **" << endl << "Press y when you're done" << endl;
    char instring[100];
    scanf("%s",instring);
  }
#endif
#endif

  // write out solution
  Yvel.Update(1.0,*Zvel,0.0);
  Ypre.Update(1.0,*Zpre,0.0);

  return 0;
}

void LINALG::SaddlePointPreconditioner::Setup(RCP<Epetra_Operator>& A,const Teuchos::ParameterList& origlist)
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


  Teuchos::RCP<Epetra_MultiVector> curvelNS = Teuchos::null;   // variables for null space
  Teuchos::RCP<Epetra_MultiVector> nextvelNS = Teuchos::null;
  Teuchos::RCP<Epetra_MultiVector> curpreNS = Teuchos::null;
  Teuchos::RCP<Epetra_MultiVector> nextpreNS = Teuchos::null;

  ///////////////// set parameter list
  RCP<Teuchos::ParameterList> spparams = Teuchos::rcp(new Teuchos::ParameterList());     // all paramaters
  RCP<Teuchos::ParameterList> velparams = Teuchos::rcp(new Teuchos::ParameterList());    // parameters (velocity specific)
  RCP<Teuchos::ParameterList> preparams = Teuchos::rcp(new Teuchos::ParameterList());    // parameters (pressure specific)

  // obtain common ML parameters from FLUID SOLVER block for coarsening from the dat file
  // we need at least "ML Parameters"."PDE equations" and "nullspace" information
  spparams->sublist("AMGBS Parameters") = params_.sublist("AMGBS Parameters"); // copy common parameters

  // first and most important: we need the number of PDE equations
  // we extract this from the Aztec Parameters and the downwind nv parameter there
  if(!params_.isSublist("Aztec Parameters")) dserror ("we expect Aztec Parameters, but there are none" );

  int nPDE = params_.sublist("Aztec Parameters").get("downwinding nv",3); // extract number of velocity dofs
  spparams->sublist("AMGBS Parameters").set("PDE equations",nPDE+1);
  spparams->sublist("AMGBS Parameters").set("null space: dimension",params_.sublist("AMGBS Parameters").get("PDE equations",3)); // copy the PDE equations as nullspace dimension

  spparams->sublist("AMGBS Parameters").set("null space: add default vectors",params_.sublist("ML Parameters").get("null space: add default vectors",false));


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
  preS_.resize(nmaxlevels_);    // smoothers
  postS_.resize(nmaxlevels_);
  Tvel_.resize(nmaxlevels_);    // transfer operators
  Tpre_.resize(nmaxlevels_);

  int nmaxcoarsedim = spparams->sublist("AMGBS Parameters").get("max coarse dimension",20);
  nVerbose_ = spparams->sublist("AMGBS Parameters").get("ML output",0);
  ndofpernode = spparams->sublist("AMGBS Parameters").get<int>("PDE equations",0);
  if(ndofpernode == 0) dserror("dof per node is zero -> error");

  nv       = ndofpernode-1;
  np       = 1;
  nlnode   = length / ndofpernode;



  /////////////////// transform Input matrix
  Ainput_ = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(A);
  if(Ainput_ != Teuchos::null)
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
    maps[0] = Teuchos::rcp(new Epetra_Map(-1,nlnode*nv,&vgid[0],0,fullmap.Comm()));
    maps[1] = Teuchos::rcp(new Epetra_Map(-1,nlnode,&pgid[0],0,fullmap.Comm()));
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
  RCP<std::vector<double> > vnewns = Teuchos::rcp(new std::vector<double>(nv*vlength,0.0));
  for (int i=0; i<nlnode; ++i)
  {
    (*vnewns)[i*nv] = 1.0;
    (*vnewns)[vlength+i*nv+1] = 1.0;
    if (nv>2) (*vnewns)[2*vlength+i*nv+2] = 1.0;
  }
  velparams->sublist("AMGBS Parameters").set("null space: vectors",&((*vnewns)[0])); // adapt default null space
  velparams->sublist("AMGBS Parameters").remove("nullspace",false);

  curvelNS = Teuchos::rcp(new Epetra_MultiVector(View,(*Ainput_)(0,0).RowMap(),&((*vnewns)[0]),(*Ainput_)(0,0).EpetraMatrix()->RowMatrixRowMap().NumMyElements(),nv));


  // pressure part: fill parameter list
  preparams->sublist("AMGBS Parameters") = spparams->sublist("AMGBS Parameters");
  preparams->sublist("AMGBS Parameters").set("PDE equations",1);               // adapt nPDE (only one pressure dof)
  preparams->sublist("AMGBS Parameters").set("null space: dimension", 1);
  const int plength = (*Ainput_)(1,1).RowMap().NumMyElements();
  RCP<std::vector<double> > pnewns = Teuchos::rcp(new std::vector<double>(plength,1.0));
  preparams->sublist("AMGBS Parameters").set("null space: vectors",&((*pnewns)[0]));
  preparams->sublist("AMGBS Parameters").remove("nullspace",false);

  curpreNS = Teuchos::rcp(new Epetra_MultiVector(View,(*Ainput_)(1,1).RowMap(),&((*pnewns)[0]),(*Ainput_)(1,1).EpetraMatrix()->RowMatrixRowMap().NumMyElements(),1));

  ////////////////// store level 0 matrices (finest level)
  int curlevel = 0;

  A11_[curlevel] = Teuchos::rcp(new SparseMatrix(Ainput_->Matrix(0,0),View));    // check me: copy or view only??
  A12_[curlevel] = Teuchos::rcp(new SparseMatrix(Ainput_->Matrix(0,1),View));
  A21_[curlevel] = Teuchos::rcp(new SparseMatrix(Ainput_->Matrix(1,0),View));
  A22_[curlevel] = Teuchos::rcp(new SparseMatrix(Ainput_->Matrix(1,1),View));

  MLAPI::Init(Teuchos::rcp(A->Comm().Clone()), true);

  for (curlevel = 0; curlevel < nmaxlevels_; ++curlevel)
  {
    /////////////////////////////////////////////////////////
    /////////////////////// AGGREGATION PROCESS
    RCP<Epetra_IntVector> velaggs = Teuchos::null; //rcp(new Epetra_IntVector(A11_[curlevel]->RowMap(),true));

    ////////////// determine aggregates using the velocity block matrix A11_[curlevel]
    RCP<AggregationMethod> aggm = AggregationMethodFactory::Create("Uncoupled",NULL);
    int naggregates_local = 0;
    int naggregates = 0;
    if(curlevel==0) velparams->sublist("AMGBS Parameters").set("Unamalgamated BlockSize",nv+1);
    else            velparams->sublist("AMGBS Parameters").set("Unamalgamated BlockSize",nv);
    velparams->sublist("AMGBS Parameters").set("Current Level",curlevel);
    naggregates = aggm->GetGlobalAggregates(A11_[curlevel]->EpetraMatrix(),velparams->sublist("AMGBS Parameters"),velaggs,naggregates_local,curvelNS);

    ////////////// transform vector with velocity aggregates to pressure block
    RCP<Epetra_IntVector> preaggs = Teuchos::rcp(new Epetra_IntVector(A22_[curlevel]->RowMap(),true));
    for(int i=0; i < preaggs->MyLength(); i++)
    {
      (*preaggs)[i] = (*velaggs)[i*nv];
    }

#ifdef WRITEOUTAGGREGATES
    // plot out aggregates
    std::stringstream fileoutstream;
    fileoutstream << "/home/wiesner/python/aggregates" << curlevel << ".vel";
    aggm->PrintIntVectorInMatlabFormat(fileoutstream.str(),*velaggs);
#endif


    /////////////////////////////////////////////////////////
    /////////////////////// CALCULATE TRANSFER OPERATORS

    ///////////// velocity transfer operators
    velparams->sublist("AMGBS Parameters").set("phase 1: max neighbour nodes", 1);
    velparams->sublist("AMGBS Parameters").set("phase 2: node attachement scheme","MaxLink");
    string velProlongSmoother = velparams->sublist("AMGBS Parameters").get("amgbs: prolongator smoother (vel)","PA-AMG");
    Tvel_[curlevel] = TransferOperatorFactory::Create(velProlongSmoother,A11_[curlevel],NULL); /* outfile */
    nextvelNS = Tvel_[curlevel]->buildTransferOperators(velaggs,naggregates_local,velparams->sublist("AMGBS Parameters"),curvelNS,0);

    //////////// pressure transfer operators
    string preProlongSmoother = preparams->sublist("AMGBS Parameters").get("amgbs: prolongator smoother (pre)","PA-AMG");
    Tpre_[curlevel] = TransferOperatorFactory::Create(preProlongSmoother,A22_[curlevel],NULL); /* outfile */
    nextpreNS = Tpre_[curlevel]->buildTransferOperators(preaggs,naggregates_local,preparams->sublist("AMGBS Parameters"),curpreNS,naggregates*nv);

    if(nVerbose_ > 4) // be verbose
    {
      cout << "Pvel[" << curlevel << "]: " << Tvel_[curlevel]->Prolongator()->EpetraMatrix()->NumGlobalRows() << " x " << Tvel_[curlevel]->Prolongator()->EpetraMatrix()->NumGlobalCols() << " (" << Tvel_[curlevel]->Prolongator()->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
      cout << "Ppre[" << curlevel << "]: " << Tpre_[curlevel]->Prolongator()->EpetraMatrix()->NumGlobalRows() << " x " << Tpre_[curlevel]->Prolongator()->EpetraMatrix()->NumGlobalCols() << " (" << Tpre_[curlevel]->Prolongator()->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;

      cout << "Rvel[" << curlevel << "]: " << Tvel_[curlevel]->Restrictor()->EpetraMatrix()->NumGlobalRows() << " x " << Tvel_[curlevel]->Restrictor()->EpetraMatrix()->NumGlobalCols() << " (" << Tvel_[curlevel]->Restrictor()->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
      cout << "Rpre[" << curlevel << "]: " << Tpre_[curlevel]->Restrictor()->EpetraMatrix()->NumGlobalRows() << " x " << Tpre_[curlevel]->Restrictor()->EpetraMatrix()->NumGlobalCols() << " (" << Tpre_[curlevel]->Restrictor()->EpetraMatrix()->NumGlobalNonzeros() << ")" << endl;
    }


    /////////////////////////// calc RAP product for next level
    A11_[curlevel+1] = Multiply(Tvel_[curlevel]->R(),*A11_[curlevel],Tvel_[curlevel]->P());
    A12_[curlevel+1] = Multiply(Tvel_[curlevel]->R(),*A12_[curlevel],Tpre_[curlevel]->P());
    A21_[curlevel+1] = Multiply(Tpre_[curlevel]->R(),*A21_[curlevel],Tvel_[curlevel]->P());
    A22_[curlevel+1] = Multiply(Tpre_[curlevel]->R(),*A22_[curlevel],Tpre_[curlevel]->P());

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
    Teuchos::ParameterList& subparams = spparams->sublist("AMGBS Parameters").sublist(stream.str());

    // copy ML Parameters or IFPACK Parameters from FLUID PRESSURE SOLVER block
    if(curlevel==0)
    {
      subparams.set("pressure correction approx: type",subparams.get("fine: type","ILU"));
      if(subparams.isSublist("IFPACK Parameters fine")) subparams.sublist("IFPACK Parameters") = subparams.sublist("IFPACK Parameters fine");
      subparams.remove("IFPACK Parameters fine",false);
      subparams.remove("IFPACK Parameters medium",false);
      subparams.remove("IFPACK Parameters coarse",false);
      subparams.remove("fine: type",false);
    }
    else
    {
      subparams.set("pressure correction approx: type",subparams.get("medium: type","ILU"));
      if(subparams.isSublist("IFPACK Parameters medium")) subparams.sublist("IFPACK Parameters") = subparams.sublist("IFPACK Parameters medium");
      subparams.remove("IFPACK Parameters fine",false);
      subparams.remove("IFPACK Parameters medium",false);
      subparams.remove("IFPACK Parameters coarse",false);
      subparams.remove("medium: type",false);
    }

    if(nVerbose_ > 8)
    {
      cout << "Braess-Sarazin smoother (level " << curlevel << ")" << endl << "parameters:" << endl << subparams << endl << endl;
    }

    preS_[curlevel]  = Teuchos::rcp(new BraessSarazin_Smoother(A11_[curlevel],A12_[curlevel],A21_[curlevel],A22_[curlevel],subparams));
    postS_[curlevel] = preS_[curlevel];//rcp(new BraessSarazin_Smoother(A11_[curlevel],A12_[curlevel],A21_[curlevel],A22_[curlevel],subparams));

    //////////////////// prepare variables for next aggregation level
    curvelNS = nextvelNS;
    curpreNS = nextpreNS;

    nlevels_ = curlevel + 1;

    //////////////////// check if aggregation is complete
    // TODO: handle aggm.getNumGlobalDirichletBlocks() in a more appropriate way
    if ((A11_[curlevel+1]->EpetraMatrix()->NumGlobalRows() + A22_[curlevel+1]->EpetraMatrix()->NumGlobalRows() - aggm->getNumGlobalDirichletBlocks()*(nv+np)) < nmaxcoarsedim)
    {
      if(nVerbose_ > 4) cout << "dim A[" << curlevel+1 << "] < " << nmaxcoarsedim << ". -> end aggregation process" << endl;
      break;
    }
  }

  //////////////////// setup coarsest smoother
  std::stringstream stream;
  stream << "braess-sarazin: list (level " << nlevels_ << ")";
  Teuchos::ParameterList& subparams = spparams->sublist("AMGBS Parameters").sublist(stream.str());

  // copy ML Parameters or IFPACK Parameters from FLUID PRESSURE SOLVER block
  subparams.set("pressure correction approx: type",subparams.get("coarse: type","ILU"));
  if(subparams.isSublist("IFPACK Parameters coarse")) subparams.sublist("IFPACK Parameters") = subparams.sublist("IFPACK Parameters coarse");
  subparams.remove("IFPACK Parameters fine",false);
  subparams.remove("IFPACK Parameters medium",false);
  subparams.remove("IFPACK Parameters coarse",false);
  subparams.remove("coarse: type",false);

  if(nVerbose_ > 8)
  {
    cout << "Braess-Sarazin smoother (level " << nlevels_ << ")" << endl << "parameters:" << endl << subparams << endl << endl;
  }

  coarsestSmoother_ = Teuchos::rcp(new BraessSarazin_Smoother(A11_[nlevels_],A12_[nlevels_],A21_[nlevels_],A22_[nlevels_],subparams));

  if(nVerbose_ > 2)
  {
    cout << "setup phase complete:" << endl;
    cout << "nlevels/maxlevels: " << nlevels_+1 << "/" << nmaxlevels_+1 << endl;
  }

  MLAPI::Finalize(true, false);


#ifdef WRITEOUTSTATISTICS
  if(outfile_)
  {
    fprintf(outfile_,"saddlepointPrecSetupTime %f\tsaddlepointPrecLevels %i\t",ttt.ElapsedTime(),nlevels_);
  }

#ifdef WRITEOUTSYMMETRY
  RCP<SparseMatrix> tmpmtx = Teuchos::rcp(new SparseMatrix(*Ainput_->Merge(),Copy));
  tmpmtx->Add(*Ainput_->Merge(),true,-1.0,1.0);
  fprintf(outfile_,"NormFrobenius %f\t",tmpmtx->NormFrobenius());
#endif
#endif


}

///////////////////////////////////////////////////////////////////
RCP<LINALG::SparseMatrix> LINALG::SaddlePointPreconditioner::Multiply(const SparseMatrix& A, const SparseMatrix& B, const SparseMatrix& C, bool bComplete)
{
  TEUCHOS_FUNC_TIME_MONITOR("SaddlePoint_Preconditioner::Multiply (with MLMultiply)");

  RCP<SparseMatrix> tmp = LINALG::MLMultiply(B,C,true);
  return LINALG::MLMultiply(A,*tmp,bComplete);
}



