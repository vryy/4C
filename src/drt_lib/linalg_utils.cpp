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
  const int ldim = (int)lm.size();
  if (ldim!=(int)lmowner.size() || ldim!=Aele.M() || ldim!=Aele.N())
    dserror("Mismatch in dimensions");

  const int myrank = A.Comm().MyPID();
  const Epetra_Map& rowmap = A.RowMap();
  const Epetra_Map& colmap = A.ColMap();

  if (A.Filled())
  {
    // loop rows of local matrix
    for (int lrow=0; lrow<ldim; ++lrow)
    {
      // check ownership of row
      if (lmowner[lrow] != myrank) continue;

      // check whether I have that global row
      int rgid = lm[lrow];
      int rlid = rowmap.LID(rgid);
#ifdef DEBUG
      if (rlid<0) dserror("Sparse matrix A does not have global row %d",rgid);
#endif

      for (int lcol=0; lcol<ldim; ++lcol)
      {
        double val = Aele(lrow,lcol);
        int cgid = lm[lcol];
        int clid = colmap.LID(cgid);
#ifdef DEBUG
        if (clid<0) dserror("Sparse matrix A does not have global column %d",cgid);
#endif
        int errone = A.SumIntoMyValues(rlid,1,&val,&clid);
        if (errone)
          dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d",errone);
      } // for (int lcol=0; lcol<ldim; ++lcol)
    } // for (int lrow=0; lrow<ldim; ++lrow)
  }
  else
  {
    // loop rows of local matrix
    for (int lrow=0; lrow<ldim; ++lrow)
    {
      // check ownership of row
      if (lmowner[lrow] != myrank) continue;

      // check whether I have that global row
      int rgid = lm[lrow];
      if (!(rowmap.MyGID(rgid))) dserror("Sparse matrix A does not have global row %d",rgid);

      for (int lcol=0; lcol<ldim; ++lcol)
      {
        double val = Aele(lrow,lcol);
        // Now that we do not rebuild the sparse mask in each step, we
        // are bound to assemble the whole thing. Zeros included.
        //if (abs(val)==0) continue; // do not assemble zeros
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
  }
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
 |  B = B*scalarB + A(transposed)*scalarA                               |
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
#define dgetrf dgetrf_
#define dgetri dgetri_
#endif
extern "C"
{
  void dsytrf(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
  void dsytri(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *info);
  void dgetrf(int *m,int *n, double *a, int *lda, int *ipiv, int* info);
  void dgetri(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
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


void LINALG::SymmetriseMatrix(Epetra_SerialDenseMatrix& A)
{
  Epetra_SerialDenseMatrix AT(A);
  AT.SetUseTranspose(true);
  // bool istranspose = AT.UseTranspose();
  A += AT;
  A.Scale(0.5);
  return;
}

/*----------------------------------------------------------------------*
 |  compute all eigenvalues and, optionally,                            |
 |  eigenvectors of a real symmetric matrix A  (public)        maf 06/07|
 *----------------------------------------------------------------------*/
void LINALG::SymmetricEigen(Epetra_SerialDenseMatrix& A,
                            Epetra_SerialDenseVector& L,
                            const int dim, const char jobz)
{
  if (A.M() != A.N()) dserror("Matrix is not square");
  if (A.M() != dim) dserror("Dimension supplied does not match matrix");
  if (L.Length() != dim) dserror("Dimension of eigenvalues does not match");

  double* a = A.A();
  double* w = L.A();
  const char uplo = {'U'};
//  char jobz = {'N'};
//  if (eigv == true) jobz = 'V';
  const int lda = A.LDA();
  const int liwork = 3+5*dim;
  vector<int> iwork(liwork);
  const int lwork = 2*dim^2 + 7*dim;
  vector<double> work(lwork);
  int info=0;

  Epetra_LAPACK lapack;

  lapack.SYEVD(jobz,uplo,dim,a,lda,w,&(work[0]),lwork,&(iwork[0]),liwork,&info);

  //SYEVD (const char JOBZ, const char UPLO, const int N, double *A, const int LDA, double *W, double *WORK, const int LWORK, int *IWORK, const int LIWORK, int *INFO) const

  if (info > 0) dserror("Lapack algorithm dsyevd failed");
  if (info < 0) dserror("Illegal value in Lapack dsyevd call");

  return;
}

/*----------------------------------------------------------------------*
 |  compute the "material tensor product" of two 2nd order tensors      |
 | (in matrix notation) and add the result to a 4th order tensor        |
 | (also in matrix notation) using the symmetry-conditions inherent to  |
 | material tensors, or tangent matrices, respectively.                 |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void LINALG::SymMatTensorMultiply(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
  // check sizes
  if (A.M() != A.N() != B.M() != B.N() != 3) dserror("2nd order tensors must be 3 by 3");
  if (C.M() != C.N() != 6) dserror("4th order tensor must be 6 by 6");

  C(0,0)= ScalarThis*C(0,0) + ScalarAB * A(0,0)*B(0,0);
  C(0,1)= ScalarThis*C(0,1) + ScalarAB * A(0,0)*B(1,1);
  C(0,2)= ScalarThis*C(0,2) + ScalarAB * A(0,0)*B(2,2);
  C(0,3)= ScalarThis*C(0,3) + ScalarAB * A(0,0)*B(1,0);
  C(0,4)= ScalarThis*C(0,4) + ScalarAB * A(0,0)*B(2,1);
  C(0,5)= ScalarThis*C(0,5) + ScalarAB * A(0,0)*B(2,0);
  
  C(1,0)= ScalarThis*C(1,0) + ScalarAB * A(1,1)*B(0,0);
  C(1,1)= ScalarThis*C(1,1) + ScalarAB * A(1,1)*B(1,1);
  C(1,2)= ScalarThis*C(1,2) + ScalarAB * A(1,1)*B(2,2);
  C(1,3)= ScalarThis*C(1,3) + ScalarAB * A(1,1)*B(1,0);
  C(1,4)= ScalarThis*C(1,4) + ScalarAB * A(1,1)*B(2,1);
  C(1,5)= ScalarThis*C(1,5) + ScalarAB * A(1,1)*B(2,0);
               
  C(2,0)= ScalarThis*C(2,0) + ScalarAB * A(2,2)*B(0,0);
  C(2,1)= ScalarThis*C(2,1) + ScalarAB * A(2,2)*B(1,1);
  C(2,2)= ScalarThis*C(2,2) + ScalarAB * A(2,2)*B(2,2);
  C(2,3)= ScalarThis*C(2,3) + ScalarAB * A(2,2)*B(1,0);
  C(2,4)= ScalarThis*C(2,4) + ScalarAB * A(2,2)*B(2,1);
  C(2,5)= ScalarThis*C(2,5) + ScalarAB * A(2,2)*B(2,0);

  C(3,0)= ScalarThis*C(3,0) + ScalarAB * A(1,0)*B(0,0);
  C(3,1)= ScalarThis*C(3,1) + ScalarAB * A(1,0)*B(1,1);
  C(3,2)= ScalarThis*C(3,2) + ScalarAB * A(1,0)*B(2,2);
  C(3,3)= ScalarThis*C(3,3) + ScalarAB * A(1,0)*B(1,0);
  C(3,4)= ScalarThis*C(3,4) + ScalarAB * A(1,0)*B(2,1);
  C(3,5)= ScalarThis*C(3,5) + ScalarAB * A(1,0)*B(2,0);
  
  C(4,0)= ScalarThis*C(4,0) + ScalarAB * A(2,1)*B(0,0);
  C(4,1)= ScalarThis*C(4,1) + ScalarAB * A(2,1)*B(1,1);
  C(4,2)= ScalarThis*C(4,2) + ScalarAB * A(2,1)*B(2,2);
  C(4,3)= ScalarThis*C(4,3) + ScalarAB * A(2,1)*B(1,0);
  C(4,4)= ScalarThis*C(4,4) + ScalarAB * A(2,1)*B(2,1);
  C(4,5)= ScalarThis*C(4,5) + ScalarAB * A(2,1)*B(2,0);
  
  C(5,0)= ScalarThis*C(5,0) + ScalarAB * A(2,0)*B(0,0);
  C(5,1)= ScalarThis*C(5,1) + ScalarAB * A(2,0)*B(1,1);
  C(5,2)= ScalarThis*C(5,2) + ScalarAB * A(2,0)*B(2,2);
  C(5,3)= ScalarThis*C(5,3) + ScalarAB * A(2,0)*B(1,0);
  C(5,4)= ScalarThis*C(5,4) + ScalarAB * A(2,0)*B(2,1);
  C(5,5)= ScalarThis*C(5,5) + ScalarAB * A(2,0)*B(2,0);

}


/*----------------------------------------------------------------------*
| invert a dense nonsymmetric matrix (public)       g.bau 03/07|
*----------------------------------------------------------------------*/
#include <Epetra_SerialDenseSolver.h>
void LINALG::NonSymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim)
{
  if (A.M() != A.N()) dserror("Matrix is not square");
  if (A.M() != dim) dserror("Dimension supplied does not match matrix");

  Epetra_SerialDenseSolver solver;
  solver.SetMatrix(A);
  int err = solver.Invert();
  if (err!=0)
    dserror("Inversion of nonsymmetric matrix failed.");

 return;
}


/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(RefCountPtr<Epetra_CrsMatrix>&   A,
                                    const RefCountPtr<Epetra_Vector> dbctoggle)
{
  RefCountPtr<Epetra_Vector> dummy = null;
  ApplyDirichlettoSystem(A,dummy,dummy,dummy,dbctoggle);
  return;
}

/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(RefCountPtr<Epetra_Vector>&      x,
                                    RefCountPtr<Epetra_Vector>&      b,
                                    const RefCountPtr<Epetra_Vector> dbcval,
                                    const RefCountPtr<Epetra_Vector> dbctoggle)
{
  RefCountPtr<Epetra_CrsMatrix> dummy = null;
  ApplyDirichlettoSystem(dummy,x,b,dbcval,dbctoggle);
  return;
}

#if 1
// The old version that modifies the sparse mask.

/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(RefCountPtr<Epetra_CrsMatrix>&   A,
                                    RefCountPtr<Epetra_Vector>&      x,
                                    RefCountPtr<Epetra_Vector>&      b,
                                    const RefCountPtr<Epetra_Vector> dbcval,
                                    const RefCountPtr<Epetra_Vector> dbctoggle)
{
  const Epetra_Vector& dbct = *dbctoggle;
  if (x != null && b != null)
  {
    Epetra_Vector&       X    = *x;
    Epetra_Vector&       B    = *b;
    const Epetra_Vector& dbcv = *dbcval;
    // set the prescribed value in x and b
    const int mylength = dbcv.MyLength();
    for (int i=0; i<mylength; ++i)
      if (dbct[i]==1.0)
      {
        X[i] = dbcv[i];
        B[i] = dbcv[i];
      }
  }

  if (A != null)
  {
    // allocate a new matrix and copy all rows that are not dirichlet
    const Epetra_Map& rowmap = A->RowMap();
    const int nummyrows     = A->NumMyRows();
    const int maxnumentries = A->MaxNumEntries();
    RefCountPtr<Epetra_CrsMatrix> Anew = LINALG::CreateMatrix(rowmap,maxnumentries);
    vector<int> indices(maxnumentries);
    vector<double> values(maxnumentries);
    for (int i=0; i<nummyrows; ++i)
    {
      int row = A->GRID(i);
      if (dbct[i]!=1.0)
      {
        int numentries;
        int err = A->ExtractGlobalRowCopy(row,maxnumentries,numentries,&values[0],&indices[0]);
#ifdef DEBUG
        if (err) dserror("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d",err);
#endif
        err = Anew->InsertGlobalValues(row,numentries,&values[0],&indices[0]);
#ifdef DEBUG
        if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
#endif
      }
      else
      {
        double one = 1.0;
        int err = Anew->InsertGlobalValues(row,1,&one,&row);
#ifdef DEBUG
        if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
#endif
      }
    }
    LINALG::Complete(*Anew);
    A = Anew;
  }

  return;
}

#else

/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                      ukue 02/07|
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(RefCountPtr<Epetra_CrsMatrix>&   A,
                                    RefCountPtr<Epetra_Vector>&      x,
                                    RefCountPtr<Epetra_Vector>&      b,
                                    const RefCountPtr<Epetra_Vector> dbcval,
                                    const RefCountPtr<Epetra_Vector> dbctoggle)
{

  const Epetra_Vector& dbct = *dbctoggle;
  if (x != null && b != null)
  {
    Epetra_Vector&       X    = *x;
    Epetra_Vector&       B    = *b;
    const Epetra_Vector& dbcv = *dbcval;
    // set the prescribed value in x and b
    const int mylength = dbcv.MyLength();
    for (int i=0; i<mylength; ++i)
      if (dbct[i]==1.0)
      {
        X[i] = dbcv[i];
        B[i] = dbcv[i];
      }
  }

  if (A != null)
  {
    const int nummyrows     = A->NumMyRows();
    for (int i=0; i<nummyrows; ++i)
    {
      //int row = A->GRID(i);
      if (dbct[i]!=1.0)
      {
        // nothing to do
      }
      else
      {
        int *indexOffset;
        int *indices;
        double *values;
        int err = A->ExtractCrsDataPointers(indexOffset, indices, values);
#ifdef DEBUG
        if (err) dserror("Epetra_CrsMatrix::ExtractCrsDataPointers returned err=%d",err);
#endif
        // zero row
        memset(&values[indexOffset[i]], 0,
               (indexOffset[i+1]-indexOffset[i])*sizeof(double));

        double one = 1.0;
        err = A->SumIntoMyValues(i,1,&one,&i);
#ifdef DEBUG
        if (err<0) dserror("Epetra_CrsMatrix::SumIntoMyValues returned err=%d",err);
#endif
      }
    }
  }

  return;
}

#endif



#endif  // #ifdef CCADISCRET
