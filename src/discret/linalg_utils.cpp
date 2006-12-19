/*!----------------------------------------------------------------------
\file linalg_utils.cpp
\brief A collection of helper methods for namespace LINALG

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "linalg_utils.H"
#include "drt_dserror.H"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "Epetra_SerialDenseSolver.h"


/*----------------------------------------------------------------------*
 |  create a Epetra_CrsMatrix  (public)                      mwgee 12/06|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_CrsMatrix> LINALG::CreateMatrix(const Epetra_Map& rowmap, const int npr)
{
  if (!rowmap.UniqueGIDs()) dserror("Row map is not unique");
  return rcp(new Epetra_CrsMatrix(Copy,rowmap,npr,false));
}
/*----------------------------------------------------------------------*
 |  create a Epetra_Vector  (public)                         mwgee 12/06|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector> LINALG::CreateVector(const Epetra_Map& rowmap, const bool init)
{
  return rcp(new Epetra_Vector(rowmap,init));
}
/*----------------------------------------------------------------------*
 |  export a Epetra_Vector  (public)                         mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::Export(const Epetra_Vector& source, Epetra_Vector& target)
{
  bool sourceunique = false;
  bool targetunique = false;
  if (source.Map().UniqueGIDs()) sourceunique = true;
  if (target.Map().UniqueGIDs()) targetunique = true;
  
  // Note:
  // source map of an import must be unique
  // target map of an export must be unique
  
  // both are unique, does not matter whether ex- or import
  if (sourceunique && targetunique)
  {
    Epetra_Export exporter(source.Map(),target.Map());
    int err = target.Export(source,exporter,Insert);
    if (err) dserror("Export using exporter returned err=%d",err);
    return;
  }
  else if (sourceunique && !targetunique)
  {
    Epetra_Import importer(target.Map(),source.Map());
    int err = target.Import(source,importer,Insert);
    if (err) dserror("Export using exporter returned err=%d",err);
    return;
  }
  else if (!sourceunique && targetunique)
  {
    Epetra_Export exporter(source.Map(),target.Map());
    int err = target.Export(source,exporter,Insert);
    if (err) dserror("Export using exporter returned err=%d",err);
    return;
  }
  else if (!sourceunique && !targetunique)
  {
    // Neither target nor source are unique - this is a problem.
    // We need a unique in between stage which we have to create artifically.
    // That's nasty.
    // As it is unclear whether this will ever be needed - do it later.
    dserror("Neither target nor source maps are unique - cannot export");
  }
  else dserror("VERY strange");
  
  return;
}

/*----------------------------------------------------------------------*
 |  assemble a matrix  (public)                              mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(Epetra_CrsMatrix& A, const Epetra_SerialDenseMatrix& Aele, 
                      const vector<int>& lm, const vector<int>& lmowner)
{
  if (A.Filled()) dserror("sparse matrix A already filled, cannot assemble");
  const int ldim = (int)lm.size();
  if (ldim!=(int)lmowner.size() || ldim!=Aele.M() || ldim!=Aele.N())
    dserror("Mismatch in dimensions");
  
  const int myrank = A.Comm().MyPID();
  
  // loop rows of local matrix
  for (int lrow=0; lrow<ldim; ++lrow)
  {
    // check ownership of row
    if (lmowner[lrow] != myrank) continue;
    
    // check whether I have that global row
    int rgid = lm[lrow];
    if (!(A.RowMap().MyGID(rgid))) dserror("Sparse matrix A does not have global row %d",rgid);
    
    for (int lcol=0; lcol<ldim; ++lcol)
    {
      double val = Aele(lrow,lcol);
      int    cgid = lm[lcol];
      int errone = A.SumIntoGlobalValues(rgid,1,&val,&cgid);
      if (errone>0)
      {
        int errtwo = A.InsertGlobalValues(rgid,1,&val,&cgid);
        if (errtwo<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code %d",errtwo);
      }
      else if (errone)
        dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d",errone);
    } // for (int lcol=0; lcol<ldim; ++lcol)
  } // for (int lrow=0; lrow<ldim; ++lrow)
  return;
}

/*----------------------------------------------------------------------*
 |  assemble a vector  (public)                              mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele, 
                const vector<int>& lm, const vector<int>& lmowner)
{
  const int ldim = (int)lm.size();
  if (ldim!=(int)lmowner.size() || ldim!=Vele.Length())
    dserror("Mismatch in dimensions");
  
  const int myrank = V.Comm().MyPID();
  
  for (int lrow=0; lrow<ldim; ++lrow)
  {
    if (lmowner[lrow] != myrank) continue;
    int rgid = lm[lrow];
    if (!V.Map().MyGID(rgid)) dserror("Sparse vector V does not have global row %d",rgid);
    int rlid = V.Map().LID(rgid);
    V[rlid] += Vele[lrow];
  } // for (int lrow=0; lrow<ldim; ++lrow)
  
  return;
}

/*----------------------------------------------------------------------*
 |  FillComplete a matrix  (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::Complete(Epetra_CrsMatrix& A)
{
  if (A.Filled()) return;
  
  int err = A.FillComplete(A.OperatorDomainMap(),A.OperatorRangeMap());
  if (err) dserror("Epetra_CrsMatrix::FillComplete(domain,range) returned err=%d",err);
  err = A.OptimizeStorage();
  if (err) dserror("Epetra_CrsMatrix::OptimizeStorage() returned err=%d",err);
  return;
}

/*----------------------------------------------------------------------*
 |  Add a sparse matrix to another  (public)                 mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::Add(const Epetra_CrsMatrix& A, 
                 const bool transposeA, 
                 const double scalarA,
                 Epetra_CrsMatrix& B, 
                 const double scalarB)
{
  if (!A.Filled()) dserror("FillComplete was not called on A");
  if (B.Filled()) dserror("FillComplete was called on B before");
  
  Epetra_CrsMatrix*               Aprime = NULL;
  EpetraExt::RowMatrix_Transpose* Atrans = NULL;
  if (transposeA)
  {
    Atrans = new EpetraExt::RowMatrix_Transpose(false,NULL,false);
    Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
  {
    Aprime = const_cast<Epetra_CrsMatrix*>(&A);
  }
  
  B.Scale(scalarB);
  
  //Loop over Aprime's rows and sum into
  int MaxNumEntries = EPETRA_MAX( Aprime->MaxNumEntries(), B.MaxNumEntries() );
  int NumEntries;
  vector<int>    Indices(MaxNumEntries);
  vector<double> Values(MaxNumEntries);

  const int NumMyRows = Aprime->NumMyRows();
  int Row, err;
  if (scalarA)
  {
    for( int i = 0; i < NumMyRows; ++i )
    {
      Row = Aprime->GRID(i);
      int ierr = Aprime->ExtractGlobalRowCopy(Row,MaxNumEntries,NumEntries,&Values[0],&Indices[0]);
      if (ierr) dserror("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d",ierr);
      if( scalarA != 1.0 )
        for( int j = 0; j < NumEntries; ++j ) Values[j] *= scalarA;
      for (int j=0; j<NumEntries; ++j)
      {
        err = B.SumIntoGlobalValues(Row,1,&Values[j],&Indices[j]);
        if (err<0 || err==2)
          err = B.InsertGlobalValues(Row,1,&Values[j],&Indices[j]);
        if (err < 0)
          dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
      }
    }
  }
  if( Atrans ) delete Atrans;
  return;
}


#ifdef LINUX_MUENCH
#define CCA_APPEND_U (1)
#endif
#ifdef CCA_APPEND_U
#define dsytrf dsytrf_
#define dsytri dsytri_
#endif
extern "C"
{
  void dsytrf(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
  void dsytri(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *info);
}

/*----------------------------------------------------------------------*
 |  invert a dense symmetric matrix  (public)                mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::SymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim)
{
  if (A.M() != A.N()) dserror("Matrix is not square");
  if (A.M() != dim) dserror("Dimension supplied does not match matrix");

  double* a = A.A();
  char uplo[5]; strncpy(uplo,"L ",2);
  vector<int> ipiv(dim);
  int lwork = 10*dim;
  vector<double> work(lwork);
  int info=0;
  int n = dim;
  int m = dim;

  dsytrf(uplo,&m,a,&n,&(ipiv[0]),&(work[0]),&lwork,&info);
  if (info) dserror("dsytrf returned info=%d",info);
  
  dsytri(uplo,&m,a,&n,&(ipiv[0]),&(work[0]),&info);
  if (info) dserror("dsytri returned info=%d",info);
  
  for (int i=0; i<dim; ++i)
    for (int j=0; j<i; ++j)
      A(j,i)=A(i,j);
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
