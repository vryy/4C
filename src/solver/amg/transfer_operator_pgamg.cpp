/*
 * transfer_operator_pgamg.cpp
 *
 *  Created on: Apr 20, 2010
 *      Author: wiesner
 */



#undef USE_SAAMG_FALLBACK

#include "transfer_operator.H"
#include "transfer_operator_saamg.H"
#include "transfer_operator_pgamg.H"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TimeMonitor.hpp"

LINALG::PGAMGTransferOperator::PGAMGTransferOperator(const RCP<SparseMatrix>& A, FILE* outfile) :
  SAAMGTransferOperator(A, outfile)
{

}

////////////////////////////////////////////////////////////
RCP<Epetra_MultiVector>  LINALG::PGAMGTransferOperator::buildTransferOperators(const RCP<Epetra_IntVector> aggs, int naggs_local, Teuchos::ParameterList& params, const RCP<Epetra_MultiVector>& ThisNS, const int domainoffset)
{
  TEUCHOS_FUNC_TIME_MONITOR("PGAMGTransferOperator::buildTransferOperators");

  ////////////// define dummy variable for next nullspace
  RCP<Epetra_MultiVector> NextNS = Teuchos::null;
  RCP<SparseMatrix> prolongator_tent = Teuchos::null;
  RCP<SparseMatrix> restrictor_tent = Teuchos::null;
  nVerbose_ = params.get("ML output",0);

  ////////////// build tentative prolongator
  GetPtent(A_->RowMap(),*aggs,naggs_local,params,*ThisNS,prolongator_tent,NextNS,domainoffset);

  ////////////// build tentative restrictor
  restrictor_tent = prolongator_tent->Transpose();

  ////////////// call internal routine for PG-AMG
  PG_AMG(A_,prolongator_tent,restrictor_tent,prolongator_,restrictor_);

  ////////////// return next level nullspace
  return NextNS;
}

////////////////////////////////////////////////////////////
void LINALG::PGAMGTransferOperator::PG_AMG(const RCP<SparseMatrix>& A, const RCP<SparseMatrix>& P_tent, const RCP<SparseMatrix>& R_tent, RCP<SparseMatrix>& P_smoothed, RCP<SparseMatrix>& R_smoothed)
{
  TEUCHOS_FUNC_TIME_MONITOR("PGAMGTransferOperator::PG_AMG");


  ////////////////// extract diagonal of A
  RCP<Epetra_Vector> diagA = Teuchos::rcp(new Epetra_Vector(A->RowMap(),true));
  A->ExtractDiagonalCopy(*diagA);
  int err = diagA->Reciprocal(*diagA);
  if(err) dserror("SaddlePointPreconditioner::PG_AMG: diagonal entries of A are 0");


  ///////////////// compute D^{-1}*A
  RCP<SparseMatrix> DinvA = Teuchos::rcp(new SparseMatrix(*A,Copy)); // ok, not the best but just works
  DinvA->LeftScale(*diagA);

  ///////////////// calculate D^{-1}*A*P0
  RCP<SparseMatrix> DinvAP0 = LINALG::MLMultiply(*DinvA,*P_tent,true);

  // TODO: drop values in DinvAP0

  // TODO: compress DinvAP0 -> DinvAP0_subset

  ///////////////// prepare variables for column-based omegas
  //int NComputedOmegas = P_tent->EpetraMatrix()->DomainMap().NumGlobalElements();  // number of col-based omegas to be computed (depends on DinvAP0_subset)

  ///////////////// compute D^{-1} * A * D^{-1} * A * P0
  RCP<SparseMatrix> DinvADinvAP0 = LINALG::MLMultiply(*DinvA,*DinvAP0,true);

  // TODO: drop small values in DinvADinvAP0

  ////////////// create map (this map makes sure, that all col based omegas can be local replicated on all processors)
#ifdef DEBUG
  if (not DinvAP0->DomainMap().UniqueGIDs())
    dserror("DinvAP0->DomainMap is not unique");
#endif
  std::vector<int> rv;
  LINALG::AllreduceEMap(rv,DinvAP0->DomainMap());
  RCP<Epetra_Map> colbasedomegamap = Teuchos::rcp(new Epetra_Map(rv.size(),rv.size(),&rv[0],0,Comm()));
#ifdef DEBUG
  if(colbasedomegamap->DistributedGlobal() == true) dserror("colbasedomegamap is distributed global?");
  if(colbasedomegamap->MinAllGID()!=DinvAP0->DomainMap().MinAllGID()) dserror("MinAllGID does not match");
  if(colbasedomegamap->MaxAllGID()!=DinvAP0->DomainMap().MaxAllGID()) dserror("MaxAllGID does not match");
  if(colbasedomegamap->NumGlobalElements()!=DinvAP0->DomainMap().NumGlobalElements()) dserror("NumGlobalElements does not match");
#endif

  ////////////// define vectors for Numerator and Denominator for calculating col based omegas
  RCP<Epetra_Vector> Numerator      = Teuchos::rcp(new Epetra_Vector(*colbasedomegamap,true));
  RCP<Epetra_Vector> Denominator    = Teuchos::rcp(new Epetra_Vector(*colbasedomegamap,true));
  RCP<Epetra_Vector> ColBasedOmegas = Teuchos::rcp(new Epetra_Vector(*colbasedomegamap,true));

  ///////////////// minimize with respect to the (D^{-1} A)' D^{-1} A norm.
  //
  //               diag( R0 (A D^{-1}' D^{-1} A' D^{-1} A' R0' )
  //  omega = ---------------------------------------------------------
  //           diag( R0 A D^{-1} A D^{-1} D^{-1} A' D^{-1} A' R0' )

  MultiplyAll(DinvAP0,DinvADinvAP0,Numerator);  // -> DinvAP0_subset
  MultiplySelfAll(DinvADinvAP0,Denominator);

#ifdef DEBUG
  int zeros_in_denominator = 0;
  for(int i=0; i<Denominator->MyLength(); i++)
  {
    if((*Denominator)[i]==0.0) zeros_in_denominator++;
  }
  if(zeros_in_denominator > 0)
  {
    cout << *Denominator << endl;
    dserror("we have %i zeros in Denominator, very suspicious",zeros_in_denominator);
  }
#endif

  ////////////// calculate column-based omega
  ColBasedOmegas->ReciprocalMultiply(1.0,*Denominator,*Numerator,0.0);

  ////////////// zero out negative column-based omegas
  int nColBasedOmegasZeros = 0;
  for(int i=0; i<ColBasedOmegas->MyLength(); i++)
  {
    if((*ColBasedOmegas)[i] < 0.0) (*ColBasedOmegas)[i] = 0.0;
    if((*ColBasedOmegas)[i] == 0.0)nColBasedOmegasZeros++;
  }

  // be verbose
  if(nVerbose_ > 6 /*&& Comm().MyPID()==0*/)
  {
    double colBasedMin = 0.0;
    double colBasedMax = 0.0;
    ColBasedOmegas->MinValue(&colBasedMin);
    ColBasedOmegas->MaxValue(&colBasedMax);
    cout << "------------------------------------------------" << endl;
    cout << "PG-AMG: damping parameters: min=" << colBasedMin << " max=" << colBasedMax << endl;
    cout << "PG-AMG: " << nColBasedOmegasZeros << " zeros out of " << ColBasedOmegas->MyLength() << " column based omegas" << endl;

  }

#ifdef USE_SAAMG_FALLBACK
  //////////////// calculate SA-AMG value for omega (for exceptions)
  double maxeig = MaxEigCG(*A,true);
  double sadampingfactor = 1.333333333333;
  double saomega = sadampingfactor/maxeig;
#endif

  //////////////// convert column based omegas to row-based omegas

  // RowBasedOmegas is a distributed vector
  RCP<Epetra_Vector> RowBasedOmegas = Teuchos::rcp(new Epetra_Vector(DinvAP0->RowMap(),true));
  RowBasedOmegas->PutScalar(-666.0);  // -666.0 -> this entry has not been set

#ifdef DEBUG
  if(!P_tent->RowMap().SameAs(DinvAP0->RowMap())) dserror("RowMaps of P_tent and DinvAP0 are not the same?");
#endif

  // some variables for statistics
  int nSpuriousZeros = 0;
  int nDirichletDofs = 0;

  // loop over local rows of DinvAP0
  for(int row=0; row<DinvAP0->EpetraMatrix()->NumMyRows(); row++)
  {
    //////////////// extract global information for local row in DinvAP0
    // we need global column ids
    int grid = DinvAP0->EpetraMatrix()->GRID(row);
    int nnz = DinvAP0->EpetraMatrix()->NumGlobalEntries(grid);
    std::vector<int> indices(nnz);
    std::vector<double> vals(nnz);
    int numEntries;
    int err = DinvAP0->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,numEntries,&vals[0],&indices[0]);
/*#ifdef DEBUG*/
    if(err!=0) dserror("Error in ExtractGlobalRowCopy");
/*#endif*/

    ///////////////// special cases
    if(nnz==0)  // empty row (should not happen, but possible because of funny MIS aggregation?)
    {
      // check if reason is an empty P_tent matrix
      if(P_tent->EpetraMatrix()->NumMyEntries(row) == 0)
      {
        // ok, this is the bad -1 aggregate thing in MIS aggregation process
#ifdef USE_SAAMG_FALLBACK
        (*RowBasedOmegas)[row] = saomega; //0.0;
#else
        (*RowBasedOmegas)[row] = 0.0;
#endif
      }
      else
      {
#ifdef DEBUG
        if(nVerbose_>=10)
        {
          if(A->EpetraMatrix()->NumMyEntries(row) == 0)
            cout << "Proc: " << Comm().MyPID() << " row(lid)=" << row << " row(gid)=" << grid << " WARNING: zero row in matrix A -> set RowOmega to zero" << endl;
          else if((*diagA)[row] == 0)
            cout << "Proc: " << Comm().MyPID() << " row(lid)=" << row << " row(gid)=" << grid << " WARNING: zero row in diagA -> set RowOmega to zero" << endl;
          else if(DinvA->EpetraMatrix()->NumMyEntries(row) == 0)
            cout << "Proc: " << Comm().MyPID() << " row(lid)=" << row << " row(gid)=" << grid << " WARNING: zero row in DinvA -> set RowOmega to zero" << endl;
          else if(DinvAP0->EpetraMatrix()->NumMyEntries(row)==0)
            cout << "Proc: " << Comm().MyPID() << " row(lid)=" << row << " row(gid)=" << grid << " WARNING: zero row in DinvAP0 -> set RowOmega to zero" << endl;
          else
            cout << "Proc: " << Comm().MyPID() << " row(lid)=" << row << " row(gid)=" << grid << " WARNING: zero row in matrix A? but why?? -> set RowOmega to zero" << endl;
          if(A->EpetraMatrix()->NumGlobalEntries(grid) == 0)
            cout << "Proc: " << Comm().MyPID() << " row(lid)=" << row << " row(gid)=" << grid << " WARNING: zero row in matrix A -> set RowOmega to zero" << endl;
          else if(DinvA->EpetraMatrix()->NumGlobalEntries(grid) == 0)
            cout << "Proc: " << Comm().MyPID() << " row(lid)=" << row << " row(gid)=" << grid << " WARNING: zero row in DinvA -> set RowOmega to zero" << endl;
        }
#endif
        nSpuriousZeros++;
#ifdef USE_SAAMG_FALLBACK
        (*RowBasedOmegas)[row] = saomega; //0.0;
#else
        (*RowBasedOmegas)[row] = 0.0;
#endif
      }
    }
    else if(nnz==1) // special treatment for dirichlet bc nodes
    {
      nDirichletDofs++;
#ifdef USE_SAAMG_FALLBACK
      (*RowBasedOmegas)[row] = saomega; //0.0;
#else
      (*RowBasedOmegas)[row] = 0.0;
#endif
    }
    else  // standard case
    {
      // find minimal ColBasedOmega within indices
      for(int j=0; j<nnz; j++)
      {
        // we need the correct "local" column index
        int col_gid = indices[j]; // this is the jth global column index in current row
        int col_omega_lid = ColBasedOmegas->Map().LID(col_gid); // this is the local id for global gid (from the ColBasedOmega view of point)
        double omega = (*ColBasedOmegas)[col_omega_lid];
        if((*RowBasedOmegas)[row]==-666.0/*DBL_MIN*/)
          (*RowBasedOmegas)[row] = omega;
        else if(omega < (*RowBasedOmegas)[row])
          (*RowBasedOmegas)[row] = omega;
      }

    }
  } // end convert col-based omegas to row-based


#ifdef DEBUG
  // check if all RowBasedOmegas are set and none entries are left
  for(int i=0; i<RowBasedOmegas->MyLength(); i++)
  {
    if((*RowBasedOmegas)[i]==-666.0)
    {
      cout << "Proc: " << Comm().MyPID() << " RowBasedOmegas[" << i << "] = -666.0 (not set)" << endl;
    }

  }
#endif

  // zero out remaining RowBasedOmegas, if they are lower than zero
  int nRowBasedOmegasZeros = 0;
  for(int i=0; i<RowBasedOmegas->MyLength(); i++)
  {
    if((*RowBasedOmegas)[i] < 0.0) (*RowBasedOmegas)[i] = 0.0;
    if((*RowBasedOmegas)[i] == 0.0) nRowBasedOmegasZeros++;
  }


  // be verbose
  if(nVerbose_ > 6 /*&& Comm().MyPID()==0*/)
  {
    double rowBasedMin = 0.0;
    double rowBasedMax = 0.0;
    RowBasedOmegas->MinValue(&rowBasedMin);
    RowBasedOmegas->MaxValue(&rowBasedMax);
    cout << "PG-AMG column->row: spurious zeros=" << nSpuriousZeros << " dirichlet bc dofs=" << nDirichletDofs << endl;
    cout << "PG-AMG: row based damping parameters: min=" << rowBasedMin << " max=" << rowBasedMax << endl;
    cout << nRowBasedOmegasZeros << " zeros out of " << RowBasedOmegas->MyLength() << " row based omegas" << endl;
    cout << "------------------------------------------------" << endl;
  }

  //////////////// compute new prolongator
  RCP<SparseMatrix> OmegaDinvA = Teuchos::rcp(new SparseMatrix(*A,Copy));
  RCP<Epetra_Vector> diagScaling = Teuchos::rcp(new Epetra_Vector(diagA->Map(),true));
  diagScaling->Multiply(1.0,*RowBasedOmegas,*diagA,0.0);
  OmegaDinvA->LeftScale(*diagScaling);
  P_smoothed = LINALG::MLMultiply(*OmegaDinvA,*P_tent,false);
  P_smoothed->Add(*P_tent,false,1.0,-1.0);
  P_smoothed->Complete(P_tent->DomainMap(),P_tent->RangeMap());

#ifdef DEBUG
  std::vector<int> smoothedzerosgids;
  int smoothedzeros = 0;
  int tentzeros = 0;
  for(int i=0; i<P_smoothed->EpetraMatrix()->NumMyRows(); i++)
  {
    int nnz = P_smoothed->EpetraMatrix()->NumMyEntries(i);
    std::vector<int> indices(nnz);
    std::vector<double> vals(nnz);
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
    std::vector<int> indices2(nnz2);
    std::vector<double> vals2(nnz2);
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

  if(nVerbose_ > 6)
    cout << "zero rows: Psmoothed=" << smoothedzeros << " Ptent=" << tentzeros << endl;

#endif

  //////////////// compute restrictor
  //R_smoothed = R_tent;
  //R_smoothed = P_smoothed->Transpose(); // geht aber wird im unsymmetrischen Teil ausgebremst
  RCP<SparseMatrix> ADinvOmega = Teuchos::rcp(new SparseMatrix(*A,Copy));
  ADinvOmega->RightScale(*diagScaling);
  R_smoothed = LINALG::MLMultiply(*R_tent,*ADinvOmega,false);
  R_smoothed->Add(*R_tent,false,1.0,-1.0);
  R_smoothed->Complete(R_tent->DomainMap(),R_tent->RangeMap());




}

void LINALG::PGAMGTransferOperator::MultiplySelfAll(const RCP<SparseMatrix>& Op,RCP<Epetra_Vector>& Column2Norm)
{
  //////////// InnerProd_local lives on current proc only, but has gids of all processors
  RCP<Epetra_Vector> Column2Norm_local = Teuchos::rcp(new Epetra_Vector(Column2Norm->Map(),true));

  //////////// collect local information
  //////////// store local scalar products in InnerProd_local vector

  // loop over local processor rows
  for(int n=0; n<Op->EpetraMatrix()->NumMyRows(); n++)
  {
    // extract values and global ids for nonzeros in current row
    int grid = Op->EpetraMatrix()->GRID(n);               // get global row id for current local row
    int nnz = Op->EpetraMatrix()->NumGlobalEntries(grid); // number of nonzero entries in current row
    std::vector<int> gindices(nnz);
    std::vector<double> gvals(nnz);
    std::vector<double> retvals(nnz);
    int gnumEntries;
    Op->EpetraMatrix()->ExtractGlobalRowCopy(grid,nnz,gnumEntries,&gvals[0],&gindices[0]);
    for(int i=0; i<(int)gvals.size(); i++)
      retvals[i] = gvals[i] * gvals[i];

    // write results back in local replicated vector with global ids
    Column2Norm_local->SumIntoGlobalValues(nnz,&retvals[0],&gindices[0]);
  }

  // safety check
  if(Column2Norm_local->GlobalLength() != Column2Norm_local->MyLength()) dserror("Global and Local length have to be the same!");

  /////////////// extract all local gids and values
  std::vector<double> localvalues(Column2Norm_local->GlobalLength(),0.0);
  std::vector<int>    globalindices(Column2Norm_local->GlobalLength(),-1);
  std::vector<double> globalvalues(Column2Norm_local->GlobalLength(),0.0);
  for(int i=0; i<(int)localvalues.size(); i++)
  {
    localvalues[i] = (*Column2Norm_local)[i];
    globalindices[i] = Column2Norm_local->Map().GID(i); // i // store indices
  }

  /////////////// sum up all local values and store it in local vector
  Comm().SumAll(&localvalues[0],&globalvalues[0],localvalues.size()); // sum up all local vectors

  /////////////// save all global information in local replicated vector Column2Norm
  Column2Norm->ReplaceGlobalValues(globalvalues.size(),&globalvalues[0],&globalindices[0]);

  //cout << "Column2Norm" << endl;
  //cout << *Column2Norm << endl;

  /////////////// free some memory
  localvalues.clear();
  globalindices.clear();
  globalvalues.clear();
}


void LINALG::PGAMGTransferOperator::MultiplyAll(const RCP<SparseMatrix>& left, const RCP<SparseMatrix>& right, RCP<Epetra_Vector>& InnerProd)
{
  if(!left->DomainMap().SameAs(right->DomainMap())) dserror("domain map of left and right does not match");
  if(!left->RowMap().SameAs(right->RowMap())) dserror("row map of left and right does not match");

  //////////// InnerProd_local lives on current proc only, but has gids of all processors
  RCP<Epetra_Vector> InnerProd_local = Teuchos::rcp(new Epetra_Vector(InnerProd->Map(),true));

  //////////// collect local information
  //////////// store local scalar products in InnerProd_local vector

  // loop over local processor rows
  for(int n=0; n<left->EpetraMatrix()->NumMyRows(); n++)
  {
    // extract values and global ids for nonzeros in current row
    int grid_left = left->EpetraMatrix()->GRID(n);               // get global row id for current local row
    int nnz_left = left->EpetraMatrix()->NumGlobalEntries(grid_left);   // number of nonzero entries in current row
    std::vector<int> gindices_left(nnz_left);
    std::vector<double> gvals_left(nnz_left);
    int gnumEntries;
    left->EpetraMatrix()->ExtractGlobalRowCopy(grid_left,nnz_left,gnumEntries,&gvals_left[0],&gindices_left[0]);

    // extract values and global ids for nonzeros in current row
    int grid_right = right->EpetraMatrix()->GRID(n);               // get global row id for current local row
    int nnz_right = right->EpetraMatrix()->NumGlobalEntries(grid_right);   // number of nonzero entries in current row
    std::vector<int> gindices_right(nnz_right);
    std::vector<double> gvals_right(nnz_right);
    right->EpetraMatrix()->ExtractGlobalRowCopy(grid_right,nnz_right,gnumEntries,&gvals_right[0],&gindices_right[0]);

    // safety check
    if(grid_left != grid_right)   dserror("grid_left and grid_right are different!");

    std::vector<double> retvals;
    std::vector<int> retindx;
    retvals.reserve(nnz_left);  // we have not more than nnz_left possible matchings of left and right gids
    retindx.reserve(nnz_left);  // we have not more than nnz_left possible matchings of left and right gids

    for(int i=0; i<(int)gvals_left.size(); i++)
    {
      for(int j=0; j<(int)gvals_right.size(); j++)
      {
        if(gindices_left[i]==gindices_right[j])
        {
          retvals.push_back(gvals_left[i]*gvals_right[j]);
          retindx.push_back(gindices_left[i]); // note: we save the gids
          break; // skip remaining gids of right operator
        }
      }
    }

    InnerProd_local->SumIntoGlobalValues((int)retvals.size(),&retvals[0],&retindx[0]); // we use the gids
  }

  // safety check
  if(InnerProd_local->GlobalLength() != InnerProd_local->MyLength()) dserror("Global and Local length have to be the same!");

  /////////////// extract all local gids and values
  std::vector<double> localvalues(InnerProd_local->GlobalLength(),0.0);
  std::vector<int>    globalindices(InnerProd_local->GlobalLength(),-1);
  std::vector<double> globalvalues(InnerProd_local->GlobalLength(),0.0);
  for(int i=0; i<(int)localvalues.size(); i++)
  {
    localvalues[i] = (*InnerProd_local)[i];
    globalindices[i] = InnerProd_local->Map().GID(i); //i; // store indices (attention: local index = global index)
  }

  /////////////// sum up all values to global
  Comm().SumAll(&localvalues[0],&globalvalues[0],localvalues.size()); // sum up all local vectors

  /////////////// save all global information in local replicated vector InnerProd
  InnerProd->ReplaceGlobalValues(globalvalues.size(),&globalvalues[0],&globalindices[0]);

  //cout << "InnerProd" << endl;
  //cout << *InnerProd << endl;

  /////////////// free some memory
  localvalues.clear();
  globalindices.clear();
  globalvalues.clear();
}


