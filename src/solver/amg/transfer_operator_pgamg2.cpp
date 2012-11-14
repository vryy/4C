/*
 * transfer_operator_pgamg2.cpp
 *
 *  Created on: Apr 23, 2010
 *      Author: wiesner
 */



#include "transfer_operator.H"
#include "transfer_operator_saamg.H"
#include "transfer_operator_pgamg2.H"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TimeMonitor.hpp"

LINALG::PGAMG2TransferOperator::PGAMG2TransferOperator(const RCP<SparseMatrix>& A, FILE* outfile) :
  SAAMGTransferOperator(A, outfile)
{

}

////////////////////////////////////////////////////////////
RCP<Epetra_MultiVector>  LINALG::PGAMG2TransferOperator::buildTransferOperators(const RCP<Epetra_IntVector> aggs, int naggs_local, Teuchos::ParameterList& params, const RCP<Epetra_MultiVector>& ThisNS, const int domainoffset)
{
  TEUCHOS_FUNC_TIME_MONITOR("PGAMG2TransferOperator::buildTransferOperators");

  ////////////// define dummy variable for next nullspace
  RCP<Epetra_MultiVector> NextNS = null;
  RCP<SparseMatrix> prolongator_tent = null;
  RCP<SparseMatrix> restrictor_tent = null;
  nVerbose_ = params.get("ML output",0);

  ////////////// build tentative prolongator
  GetPtent(A_->RowMap(),*aggs,naggs_local,params,*ThisNS,prolongator_tent,NextNS,domainoffset);

  ////////////// build tentative restrictor
  restrictor_tent = prolongator_tent->Transpose();

  /////////////// extract diagonal of A
  RCP<Epetra_Vector> diagA = Teuchos::rcp(new Epetra_Vector(A_->RowMap(),true));
  A_->ExtractDiagonalCopy(*diagA);
  int err = diagA->Reciprocal(*diagA);
  if(err) dserror("PGAMGTransferOperator::buildTransferOperators: diagonal entries of A are 0");

  ///////////////// transpose A
  RCP<SparseMatrix> At = A_->Transpose();

  ////////////// compute AP0 and AtP0
  RCP<SparseMatrix> AP0 = LINALG::MLMultiply(*A_,*prolongator_tent,true);
  RCP<SparseMatrix> AtP0 = LINALG::MLMultiply(*At,*prolongator_tent,true);

  ////////////// calculate col based omegas for P and row based omegas for R
  RCP<Epetra_Vector> ColBasedOmegaP = CalculateColBasedOmegaP(A_,diagA,prolongator_tent,AP0);
  RCP<Epetra_Vector> RowBasedOmegaR = CalculateRowBasedOmegaR(At,diagA,prolongator_tent,AtP0);

/*#ifdef DEBUG
  RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(ColBasedOmegaP->Map(),true));
  temp->Update(1.0,*ColBasedOmegaP,0.0);
  temp->Update(-1.0,*RowBasedOmegaR,1.0);
  cout << "difference of omegas" << endl;
  cout << *temp << endl;
#endif*/

  ////////////// calculate row based omegas for P and col based omegas for R
  RCP<Epetra_Vector> RowBasedOmegaP = TransformColBased2RowBasedOmegas(AP0,ColBasedOmegaP,0.0);
  RCP<Epetra_Vector> ColBasedOmegaR = TransformColBased2RowBasedOmegas(AtP0,RowBasedOmegaR,0.0);

/*#ifdef DEBUG
  RCP<Epetra_Vector> temp2 = Teuchos::rcp(new Epetra_Vector(RowBasedOmegaP->Map(),true));
  temp2->Update(1.0,*RowBasedOmegaP,0.0);
  temp2->Update(-1.0,*ColBasedOmegaR,1.0);
  cout << "difference of omegas" << endl;
  cout << *temp2 << endl;
#endif*/

  ////////////// calculate prolongator
  RCP<SparseMatrix> OmegaDinvA = Teuchos::rcp(new SparseMatrix(*A_,Copy));
  RCP<Epetra_Vector> diagScaling = Teuchos::rcp(new Epetra_Vector(diagA->Map(),true));
  diagScaling->Multiply(1.0,*RowBasedOmegaP,*diagA,0.0);
  OmegaDinvA->LeftScale(*diagScaling);
  prolongator_ = LINALG::MLMultiply(*OmegaDinvA,*prolongator_tent,false);
  prolongator_->Add(*prolongator_tent,false,1.0,-1.0);
  prolongator_->Complete(prolongator_tent->DomainMap(),prolongator_tent->RangeMap());

  /////////////// calculate restrictor
  RCP<SparseMatrix> ADinvOmega = Teuchos::rcp(new SparseMatrix(*A_,Copy));
  diagScaling->Multiply(1.0,*diagA,*ColBasedOmegaR,0.0);
  ADinvOmega->RightScale(*diagScaling);
  restrictor_ = LINALG::MLMultiply(*restrictor_tent,*ADinvOmega,false);
  restrictor_->Add(*restrictor_tent,false,1.0,-1.0);
  restrictor_->Complete(restrictor_tent->DomainMap(),restrictor_tent->RangeMap());

  ////////////// return next level nullspace
  return NextNS;
}


RCP<Epetra_Vector> LINALG::PGAMG2TransferOperator::CalculateColBasedOmegaP(const RCP<SparseMatrix>& A,const RCP<Epetra_Vector>& Dinv,const RCP<SparseMatrix>& Ptent, const RCP<SparseMatrix>& AP0)
{
  // return value: colbased omegaP vector

  /////////////// create map for col based omegas (that is replicated on all procs)
  std::vector<int> rv;
  LINALG::AllreduceEMap(rv,Ptent->DomainMap());
  RCP<Epetra_Map> colbasedomegamap = Teuchos::rcp(new Epetra_Map(rv.size(),rv.size(),&rv[0],0,Comm()));

#ifdef DEBUG
  if(not Ptent->DomainMap().UniqueGIDs()) dserror("Domain map of Ptent is not unique");
  if(colbasedomegamap->DistributedGlobal() == true) dserror("colbasedomegamap ist distributed global and not locally replicated");
  if(colbasedomegamap->MinAllGID()!=Ptent->DomainMap().MinAllGID()) dserror("MinAllGID does not match");
  if(colbasedomegamap->MaxAllGID()!=Ptent->DomainMap().MaxAllGID()) dserror("MaxAllGID does not match");
  if(colbasedomegamap->NumGlobalElements()!=Ptent->DomainMap().NumGlobalElements()) dserror("NumGlobalElements does not match");
#endif


  /////////////// compute ADinv
  RCP<SparseMatrix> ADinv = Teuchos::rcp(new SparseMatrix(*A,Copy));
  ADinv->RightScale(*Dinv);

  /////////////// calculate ADinvA
  RCP<SparseMatrix> ADinvA = LINALG::MLMultiply(*ADinv,*A,true);

  /////////////// calculate ADinvAPtent
  RCP<SparseMatrix> ADinvAP0 = LINALG::MLMultiply(*ADinvA,*Ptent,true);

  /////////////// define vectors using colbasedomegamap
  RCP<Epetra_Vector> Numerator = Teuchos::rcp(new Epetra_Vector(*colbasedomegamap,true));
  RCP<Epetra_Vector> Denominator = Teuchos::rcp(new Epetra_Vector(*colbasedomegamap,true));
  RCP<Epetra_Vector> ColBasedOmegasP = Teuchos::rcp(new Epetra_Vector(*colbasedomegamap,true));

  /////////////// minimize with respect to A' A norm
  //
  //           diag( P0' A' A Dinv A P0 )
  // omega = ------------------------------------
  //          diag( P0' A' Dinv A' A Dinv A P0 )

  MultiplyAll(AP0,ADinvAP0,Numerator);
  MultiplySelfAll(ADinvAP0,Denominator);

#ifdef DEBUG
  int zeros_in_denominator = 0;
  for(int i=0; i<Denominator->MyLength(); i++)
  {
    if((*Denominator)[i]==0.0) zeros_in_denominator++;
  }
  if(zeros_in_denominator > 0)
  {
    cout << *Denominator << endl;
    dserror("we have %i zeros in Denominator, very suspicious", zeros_in_denominator);
  }
#endif

  ////////////// calculate column-based omega
  ColBasedOmegasP->ReciprocalMultiply(1.0,*Denominator,*Numerator,0.0);

  ////////////// zero out negative column-based omegas
  int nColBasedOmegasZeros = 0;
  for(int i=0; i<ColBasedOmegasP->MyLength(); i++)
  {
    if((*ColBasedOmegasP)[i] < 0.0) (*ColBasedOmegasP)[i] = 0.0;
    if((*ColBasedOmegasP)[i] == 0.0) nColBasedOmegasZeros++;
  }

  // be verbose
  if(nVerbose_ > 6)
  {
    double colBasedMin = 0.0;
    double colBasedMax = 0.0;
    ColBasedOmegasP->MinValue(&colBasedMin);
    ColBasedOmegasP->MaxValue(&colBasedMax);
    cout << "PG-AMG: Proc: " << Comm().MyPID() << " min=" << colBasedMin << " max=" << colBasedMax << endl;
    cout << "PG-AMG: Proc: " << Comm().MyPID() << " " << nColBasedOmegasZeros << " zeros out of " << ColBasedOmegasP->MyLength() << " column based omegas" << endl;
  }

  return ColBasedOmegasP;
}

RCP<Epetra_Vector> LINALG::PGAMG2TransferOperator::CalculateRowBasedOmegaR(const RCP<SparseMatrix>& At,const RCP<Epetra_Vector>& Dinv,const RCP<SparseMatrix>& Ptent,const RCP<SparseMatrix>& AtP0)
{

  /////////////// compute AtDinv
  RCP<SparseMatrix> AtDinv = Teuchos::rcp(new SparseMatrix(*At,Copy));
  AtDinv->RightScale(*Dinv);

  /////////////// calculate AtDinvAt
  RCP<SparseMatrix> AtDinvAt = LINALG::MLMultiply(*AtDinv,*At,true);

  /////////////// calculate AtDinvAtPtent
  RCP<SparseMatrix> AtDinvAtP0 = LINALG::MLMultiply(*AtDinvAt,*Ptent,true);

  /////////////// create map for col based omegas (that is replicated on all procs)
  std::vector<int> rv;
  LINALG::AllreduceEMap(rv,Ptent->DomainMap());
  RCP<Epetra_Map> rowbasedomegamar = Teuchos::rcp(new Epetra_Map(rv.size(),rv.size(),&rv[0],0,Comm()));

#ifdef DEBUG
  if(not Ptent->DomainMap().UniqueGIDs()) dserror("Domain map of Ptent is not unique");
  if(rowbasedomegamar->DistributedGlobal() == true) dserror("colbasedomegamap ist distributed global and not locally replicated");
  if(rowbasedomegamar->MinAllGID()!=Ptent->DomainMap().MinAllGID()) dserror("MinAllGID does not match");
  if(rowbasedomegamar->MaxAllGID()!=Ptent->DomainMap().MaxAllGID()) dserror("MaxAllGID does not match");
  if(rowbasedomegamar->NumGlobalElements()!=Ptent->DomainMap().NumGlobalElements()) dserror("NumGlobalElements does not match");
#endif

  /////////////// define vectors using colbasedomegamap
  RCP<Epetra_Vector> Numerator = Teuchos::rcp(new Epetra_Vector(*rowbasedomegamar,true));
  RCP<Epetra_Vector> Denominator = Teuchos::rcp(new Epetra_Vector(*rowbasedomegamar,true));
  RCP<Epetra_Vector> RowBasedOmegasR = Teuchos::rcp(new Epetra_Vector(*rowbasedomegamar,true));

  /////////////// minimize with respect to A' A norm
  //
  //           diag( P0' A A' Dinv A' P0 )
  // omega = ------------------------------------
  //          diag( P0' A Dinv A A' Dinv A' P0 )

  MultiplyAll(AtP0,AtDinvAtP0,Numerator);
  MultiplySelfAll(AtDinvAtP0,Denominator);

#ifdef DEBUG
  int zeros_in_denominator = 0;
  for(int i=0; i<Denominator->MyLength(); i++)
  {
    if((*Denominator)[i]==0.0) zeros_in_denominator++;
  }
  if(zeros_in_denominator > 0)
  {
    cout << *Denominator << endl;
    dserror("we have %i zeros in Denominator, very suspicious", zeros_in_denominator);
  }
#endif

  ////////////// calculate row-based omega for restrictor
  RowBasedOmegasR->ReciprocalMultiply(1.0,*Denominator,*Numerator,0.0);

  ////////////// zero out negative column-based omegas
  int nRowBasedOmegasZeros = 0;
  for(int i=0; i<RowBasedOmegasR->MyLength(); i++)
  {
    if((*RowBasedOmegasR)[i] < 0.0) (*RowBasedOmegasR)[i] = 0.0;
    if((*RowBasedOmegasR)[i] == 0.0) nRowBasedOmegasZeros++;
  }

  // be verbose
  if(nVerbose_ > 6)
  {
    double rowBasedMin = 0.0;
    double rowBasedMax = 0.0;
    RowBasedOmegasR->MinValue(&rowBasedMin);
    RowBasedOmegasR->MaxValue(&rowBasedMax);
    cout << "PG-AMG: Proc: " << Comm().MyPID() << " min=" << rowBasedMin << " max=" << rowBasedMax << endl;
    cout << "PG-AMG: Proc: " << Comm().MyPID() << " " << nRowBasedOmegasZeros << " zeros out of " << RowBasedOmegasR->MyLength() << " row based restrictor omegas" << endl;
  }

  return RowBasedOmegasR;
}

RCP<Epetra_Vector> LINALG::PGAMG2TransferOperator::TransformColBased2RowBasedOmegas(const RCP<SparseMatrix>& AP, const RCP<Epetra_Vector> colbasedomegas, double fallbackvalue)
{
  /////////////// define row based omega vector (distributed over row map of AP matrix)
  RCP<Epetra_Vector> RowBasedOmegas = Teuchos::rcp(new Epetra_Vector(AP->RowMap(),true));

  /////////////// mark all row based omegas as not defined (that is: set all entries to -666.0)
  RowBasedOmegas->PutScalar(-666.0);

  /////////////// define some variables for statistics
  int nDirichletDofs = 0;
  int nSpuriousZeros = 0;

  // loop over all local rows of AP
  for(int row=0; row<AP->EpetraMatrix()->NumMyRows(); row++)
  {
    /////////////// extract global information for current local row in AP
    int grid = AP->EpetraMatrix()->GRID(row); // global row id
    int nnz  = AP->EpetraMatrix()->NumGlobalEntries(grid);  // number of nonzeros in current row

    /////////////// special cases
    if(nnz==0)
    { // this should not be happen (but is possible because of the bad MIS aggregation)
      nSpuriousZeros++;
      (*RowBasedOmegas)[row] = fallbackvalue;
    }
    else if(nnz==1)
    {
      // dirichlet bc nodes
      nDirichletDofs++;
      (*RowBasedOmegas)[row] = fallbackvalue;
    }
    else
    {
      ////////////// standard case
      std::vector<int> indices(nnz);
      std::vector<double> values(nnz);
      int numEntries;
      int err = AP->EpetraMatrix()->ExtractGlobalRowCopy(grid,
                                                         nnz,
                                                         numEntries,
                                                         &values[0],
                                                         &indices[0]);
//#ifdef DEBUG
      if(err!=0) dserror("Error in ExtractGlobalRowCopy");
//#endif

      ///////////// find corresponding minimal colbased omega within indices
      for(int j=0; j<nnz; j++)
      {
        // we  need the correct "local" global column index
        int col_gid = indices[j]; // this is the jth global column index in the current row
        int col_omega_lid = colbasedomegas->Map().LID(col_gid);// this is the local id for the global gid (from the colbasedomega view of point)
        double omega = (*colbasedomegas)[col_omega_lid];
        if((*RowBasedOmegas)[row] == -666.0 || omega < (*RowBasedOmegas)[row])
          (*RowBasedOmegas)[row] = omega;
      } // end loop over all columns

    }

  } // end loop over all rows of AP

#ifdef DEBUG
  // check if all row based omegas are set and no entries are left
  for(int i=0; i<RowBasedOmegas->MyLength(); i++)
  {
    if((*RowBasedOmegas)[i]==-666.0) // this should not happen
    {
      cout << "Proc: " << Comm().MyPID() << " RowBasedOmegas[" << i<< "] = -666.0 (not set)" << endl;
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
  if(nVerbose_ > 6)
  {
    double rowBasedMin = 0.0;
    double rowBasedMax = 0.0;
    RowBasedOmegas->MinValue(&rowBasedMin);
    RowBasedOmegas->MaxValue(&rowBasedMax);
    cout << "PG-AMG: Proc: " << Comm().MyPID() << " min=" << rowBasedMin << " max=" << rowBasedMax << endl;
    cout << "PG-AMG: Proc: " << Comm().MyPID() << " " << nRowBasedOmegasZeros << " zeros out of " << RowBasedOmegas->MyLength() << " row based omegas" << endl;
  }

  return RowBasedOmegas;
}

void LINALG::PGAMG2TransferOperator::MultiplySelfAll(const RCP<SparseMatrix>& Op,RCP<Epetra_Vector>& Column2Norm)
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


void LINALG::PGAMG2TransferOperator::MultiplyAll(const RCP<SparseMatrix>& left, const RCP<SparseMatrix>& right, RCP<Epetra_Vector>& InnerProd)
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


