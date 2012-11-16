/*
 * transfer_operator_tentative.cpp
 *
 *  Created on: Apr 20, 2010
 *      Author: wiesner
 */


#include "transfer_operator.H"
#include "transfer_operator_tentative.H"

LINALG::TentativeTransferOperator::TentativeTransferOperator(const RCP<SparseMatrix>& A, FILE* outfile) :
  TransferOperator(A,outfile)
{

}

RCP<Epetra_MultiVector>  LINALG::TentativeTransferOperator::buildTransferOperators(const RCP<Epetra_IntVector> aggs, int naggs_local, Teuchos::ParameterList& params, const RCP<Epetra_MultiVector>& ThisNS, const int domainoffset)
{
  ////////////// define dummy variable for next nullspace
  RCP<Epetra_MultiVector> NextNS = Teuchos::null;

  ////////////// build tentative prolongator
  GetPtent(A_->RowMap(),*aggs,naggs_local,params,*ThisNS,prolongator_,NextNS,domainoffset);

  ////////////// build tentative restrictor
  restrictor_ = prolongator_->Transpose();

  ////////////// return next level nullspace
  return NextNS;
}

void LINALG::TentativeTransferOperator::GetPtent(const Epetra_Map& rowmap, const Epetra_IntVector aggvec, int naggs, ParameterList& List, const Epetra_MultiVector& ThisNS, RCP<SparseMatrix>& Ptent, RCP<Epetra_MultiVector>& NextNS, const int domainoffset)
{
  const int nsdim = List.get<int>("null space: dimension",-1);
  if(nsdim <= 0) dserror("null space dimension not given");
  const int mylength = rowmap.NumMyElements();

  //////////////// build a domain map for Ptent
#if 0
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
#else
  // we suppose the aggids to be ordered sequential on each proc
  // search aggid with minimal value
  int firstagg = 100000000;    // agg with minimal aggid on current proc
  int offset = -1;
  for (int i=0; i<mylength; ++i)
  {
    if (aggvec[i]>=0 && aggvec[i]<firstagg)
    {
      firstagg = aggvec[i];
    }
  }
  offset = firstagg * nsdim;
#endif

  // generate gids for coarse grid
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
    aggdofs.insert(std::pair<int,std::vector<int> >(firstagg+i,gids));  // is meant to contain all dof gids for an aggregate
  }
  for(int i=0; i<mylength; ++i)
  {
    if(aggvec[i] < 0) continue;   // this agg doesn't belong to current proc
    vector<int>& gids = aggdofs[aggvec[i]];
    gids.push_back(aggvec.Map().GID(i));  // add dof gid to gids for current aggregate
  }

  //////////////// coarse level nullspace to be filled
  NextNS = Teuchos::rcp(new Epetra_MultiVector(pdomainmap,nsdim,true));
  Epetra_MultiVector& nextns = * NextNS;

  //////////////// create Ptent
  Ptent = Teuchos::rcp(new SparseMatrix(rowmap, nsdim));   // create new Ptent matrix

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

    std::vector<double> work(lwork);
    std::vector<double> tau(k);
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




