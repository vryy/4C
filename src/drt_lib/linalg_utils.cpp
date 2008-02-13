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

#include <algorithm>
#include <numeric>
#include <vector>

#include "linalg_utils.H"
#include "linalg_systemmatrix.H"
#include "drt_dserror.H"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_RowMatrixTransposer.h"


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
void LINALG::Export(const Epetra_MultiVector& source, Epetra_MultiVector& target)
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
 |  assemble a matrix  (public)                               popp 01/08|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(Epetra_CrsMatrix& A, const Epetra_SerialDenseMatrix& Aele,
                      const vector<int>& lmrow, const vector<int>& lmrowowner,
                      const vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();
  if (lrowdim!=(int)lmrowowner.size() || lrowdim!=Aele.M() || lcoldim!=Aele.N())
    dserror("Mismatch in dimensions");

  const int myrank = A.Comm().MyPID();
  const Epetra_Map& rowmap = A.RowMap();

  // this 'Assemble' is not implemented for a Filled() matrix A
  if (A.Filled()) dserror("Sparse matrix A is already Filled()");

  else
  {
    // loop rows of local matrix
    for (int lrow=0; lrow<lrowdim; ++lrow)
    {
      // check ownership of row
      if (lmrowowner[lrow] != myrank) continue;

      // check whether I have that global row
      int rgid = lmrow[lrow];
      if (!(rowmap.MyGID(rgid))) dserror("Sparse matrix A does not have global row %d",rgid);

      for (int lcol=0; lcol<lcoldim; ++lcol)
      {
        double val = Aele(lrow,lcol);
        int cgid = lmcol[lcol];

        // Now that we do not rebuild the sparse mask in each step, we
        // are bound to assemble the whole thing. Zeros included.
        int errone = A.SumIntoGlobalValues(rgid,1,&val,&cgid);
        if (errone>0)
        {
          int errtwo = A.InsertGlobalValues(rgid,1,&val,&cgid);
          if (errtwo<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code %d",errtwo);
        }
        else if (errone)
          dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d",errone);
      } // for (int lcol=0; lcol<lcoldim; ++lcol)
    } // for (int lrow=0; lrow<lrowdim; ++lrow)
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
 |  assemble a vector into MultiVector (public)              mwgee 01/08|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
                const vector<int>& lm, const vector<int>& lmowner)
{
  LINALG::Assemble(*(V(n)),Vele,lm,lmowner);
  return;
}

/*----------------------------------------------------------------------*
 |  FillComplete a matrix  (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::Complete(Epetra_CrsMatrix& A)
{
  if (A.Filled()) return;
  int err = A.FillComplete(A.OperatorDomainMap(),A.OperatorRangeMap(),true);
  if (err) dserror("Epetra_CrsMatrix::FillComplete(domain,range) returned err=%d",err);
  return;
}

/*----------------------------------------------------------------------*
 |  FillComplete a matrix  (public)                          mwgee 01/08|
 *----------------------------------------------------------------------*/
void  LINALG::Complete(Epetra_CrsMatrix& A, const Epetra_Map& domainmap, const Epetra_Map& rangemap)
{
  if (A.Filled()) return;
  int err = A.FillComplete(domainmap,rangemap,true);
  if (err) dserror("Epetra_CrsMatrix::FillComplete(domain,range) returned err=%d",err);
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

  Epetra_CrsMatrix* Aprime = NULL;
  RCP<EpetraExt::RowMatrix_Transpose> Atrans = null;
  if (transposeA)
  {
    Atrans = rcp(new EpetraExt::RowMatrix_Transpose(false,NULL,false));
    Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
  {
    Aprime = const_cast<Epetra_CrsMatrix*>(&A);
  }

  if (scalarB != 1.0) B.Scale(scalarB);
  if (scalarB == 0.0) B.PutScalar(0.0);

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
  return;
}

/*----------------------------------------------------------------------*
 | Transpose matrix A                                         popp 02/08|
 *----------------------------------------------------------------------*/
RCP<Epetra_CrsMatrix> LINALG::Transpose(const Epetra_CrsMatrix& A)
{
  if (!A.Filled()) dserror("FillComplete was not called on A");

  if (!A.Filled()) dserror("FillComplete was not called on A");

  RCP<EpetraExt::RowMatrix_Transpose> Atrans =
  		rcp(new EpetraExt::RowMatrix_Transpose(false,NULL,false));
  Epetra_CrsMatrix* Aprime =
  		&(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));

  
  return rcp(new Epetra_CrsMatrix(*Aprime));
}

/*----------------------------------------------------------------------*
 | Multiply matrices A*B                                     mwgee 01/06|
 *----------------------------------------------------------------------*/
RCP<Epetra_CrsMatrix> LINALG::Multiply(const Epetra_CrsMatrix& A, bool transA,
                                       const Epetra_CrsMatrix& B, bool transB)
{
  // make sure FillComplete was called on the matrices
  if (!A.Filled()) dserror("A has to be FillComplete");
  if (!B.Filled()) dserror("B has to be FillComplete");

  // create resultmatrix with correct rowmap
  Epetra_CrsMatrix* C = NULL;
  if (!transA)
    C = new Epetra_CrsMatrix(Copy,A.OperatorRangeMap(),20,false);
  else
    C = new Epetra_CrsMatrix(Copy,A.OperatorDomainMap(),20,false);

  int err = EpetraExt::MatrixMatrix::Multiply(A,transA,B,transB,*C);
  if (err) dserror("EpetraExt::MatrixMatrix::Multiply returned err = &d",err);

  return rcp(C);
}

/*----------------------------------------------------------------------*
 | Multiply matrices A*B                                     mwgee 02/08|
 *----------------------------------------------------------------------*/
RCP<Epetra_CrsMatrix> LINALG::Multiply(const Epetra_CrsMatrix& A, bool transA,
                                       const Epetra_CrsMatrix& B, bool transB,
                                       const Epetra_CrsMatrix& C, bool transC)
{
  RCP<Epetra_CrsMatrix> tmp = LINALG::Multiply(B,transB,C,transC);
  return LINALG::Multiply(A,transA,*tmp,false);
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


/*----------------------------------------------------------------------*
 |                                             (public)        maf 06/07|
 *----------------------------------------------------------------------*/
void LINALG::SymmetriseMatrix(Epetra_SerialDenseMatrix& A)
{
#if 1
  const int n = A.N();
  if (n != A.M()) dserror("Cannot symmetrize non-square matrix");
  // do not make deep copy of A, matrix addition and full scaling just to sym it
  // gee
  for (int i=0; i<n; ++i)
    for (int j=i+1; j<n; ++j)
    {
      const double aver = 0.5*(A(i,j)+A(j,i));
      A(i,j) = A(j,i) = aver;
    }
#else
  Epetra_SerialDenseMatrix AT(A);
  AT.SetUseTranspose(true);
  // bool istranspose = AT.UseTranspose();
  A += AT;
  A.Scale(0.5);
#endif
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
 |  compute the "material tensor product" A x B of two 2nd order tensors|
 | (in matrix notation) and add the result to a 4th order tensor        |
 | (also in matrix notation) using the symmetry-conditions inherent to  |
 | material tensors, or tangent matrices, respectively.                 |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void LINALG::ElastSymTensorMultiply(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");

  // everything in Voigt-Notation
  Epetra_SerialDenseMatrix AVoigt(6,1);
  Epetra_SerialDenseMatrix BVoigt(6,1);
  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  // Voigts vector notation implies 2 times ()12 ()23 ()13 !
  //AVoigt(3,0) = 2.0*A(1,0); AVoigt(4,0) = 2.0*A(2,1); AVoigt(5,0) = 2.0*A(2,0);
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);
  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  // Voigts vector notation implies 2 times ()12 ()23 ()13 !
  //BVoigt(3,0) = 2.0*B(1,0); BVoigt(4,0) = 2.0*B(2,1); BVoigt(5,0) = 2.0*B(2,0);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);
  C.Multiply('N','T',ScalarAB,AVoigt,BVoigt,ScalarThis);

  // this is explicitly what the former .Multiply does:
//  C(0,0)= ScalarThis*C(0,0) + ScalarAB * A(0,0)*B(0,0);
//  C(0,1)= ScalarThis*C(0,1) + ScalarAB * A(0,0)*B(1,1);
//  C(0,2)= ScalarThis*C(0,2) + ScalarAB * A(0,0)*B(2,2);
//  C(0,3)= ScalarThis*C(0,3) + ScalarAB * A(0,0)*B(1,0);
//  C(0,4)= ScalarThis*C(0,4) + ScalarAB * A(0,0)*B(2,1);
//  C(0,5)= ScalarThis*C(0,5) + ScalarAB * A(0,0)*B(2,0);
//
//  C(1,0)= ScalarThis*C(1,0) + ScalarAB * A(1,1)*B(0,0);
//  C(1,1)= ScalarThis*C(1,1) + ScalarAB * A(1,1)*B(1,1);
//  C(1,2)= ScalarThis*C(1,2) + ScalarAB * A(1,1)*B(2,2);
//  C(1,3)= ScalarThis*C(1,3) + ScalarAB * A(1,1)*B(1,0);
//  C(1,4)= ScalarThis*C(1,4) + ScalarAB * A(1,1)*B(2,1);
//  C(1,5)= ScalarThis*C(1,5) + ScalarAB * A(1,1)*B(2,0);
//
//  C(2,0)= ScalarThis*C(2,0) + ScalarAB * A(2,2)*B(0,0);
//  C(2,1)= ScalarThis*C(2,1) + ScalarAB * A(2,2)*B(1,1);
//  C(2,2)= ScalarThis*C(2,2) + ScalarAB * A(2,2)*B(2,2);
//  C(2,3)= ScalarThis*C(2,3) + ScalarAB * A(2,2)*B(1,0);
//  C(2,4)= ScalarThis*C(2,4) + ScalarAB * A(2,2)*B(2,1);
//  C(2,5)= ScalarThis*C(2,5) + ScalarAB * A(2,2)*B(2,0);
//
//  C(3,0)= ScalarThis*C(3,0) + ScalarAB * A(1,0)*B(0,0);
//  C(3,1)= ScalarThis*C(3,1) + ScalarAB * A(1,0)*B(1,1);
//  C(3,2)= ScalarThis*C(3,2) + ScalarAB * A(1,0)*B(2,2);
//  C(3,3)= ScalarThis*C(3,3) + ScalarAB * A(1,0)*B(1,0);
//  C(3,4)= ScalarThis*C(3,4) + ScalarAB * A(1,0)*B(2,1);
//  C(3,5)= ScalarThis*C(3,5) + ScalarAB * A(1,0)*B(2,0);
//
//  C(4,0)= ScalarThis*C(4,0) + ScalarAB * A(2,1)*B(0,0);
//  C(4,1)= ScalarThis*C(4,1) + ScalarAB * A(2,1)*B(1,1);
//  C(4,2)= ScalarThis*C(4,2) + ScalarAB * A(2,1)*B(2,2);
//  C(4,3)= ScalarThis*C(4,3) + ScalarAB * A(2,1)*B(1,0);
//  C(4,4)= ScalarThis*C(4,4) + ScalarAB * A(2,1)*B(2,1);
//  C(4,5)= ScalarThis*C(4,5) + ScalarAB * A(2,1)*B(2,0);
//
//  C(5,0)= ScalarThis*C(5,0) + ScalarAB * A(2,0)*B(0,0);
//  C(5,1)= ScalarThis*C(5,1) + ScalarAB * A(2,0)*B(1,1);
//  C(5,2)= ScalarThis*C(5,2) + ScalarAB * A(2,0)*B(2,2);
//  C(5,3)= ScalarThis*C(5,3) + ScalarAB * A(2,0)*B(1,0);
//  C(5,4)= ScalarThis*C(5,4) + ScalarAB * A(2,0)*B(2,1);
//  C(5,5)= ScalarThis*C(5,5) + ScalarAB * A(2,0)*B(2,0);

  return;
}

/*----------------------------------------------------------------------*
 | compute the "material tensor product" (A x B + B x A) of two 2nd     |
 | order tensors                                                        |
 | (in matrix notation) and add the result to a 4th order tensor        |
 | (also in matrix notation) using the symmetry-conditions inherent to  |
 | material tensors, or tangent matrices, respectively.                 |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void LINALG::ElastSymTensorMultiplyAddSym(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");

  // everything in Voigt-Notation
  Epetra_SerialDenseMatrix AVoigt(6,1);
  Epetra_SerialDenseMatrix BVoigt(6,1);
  AVoigt(0,0) = A(0,0); AVoigt(1,0) = A(1,1); AVoigt(2,0) = A(2,2);
  // Voigts vector notation implies 2 times ()12 ()23 ()13 !
  //AVoigt(3,0) = 2.0*A(1,0); AVoigt(4,0) = 2.0*A(2,1); AVoigt(5,0) = 2.0*A(2,0);
  AVoigt(3,0) = A(1,0); AVoigt(4,0) = A(2,1); AVoigt(5,0) = A(2,0);
  BVoigt(0,0) = B(0,0); BVoigt(1,0) = B(1,1); BVoigt(2,0) = B(2,2);
  // Voigts vector notation implies 2 times ()12 ()23 ()13 !
  //BVoigt(3,0) = 2.0*B(1,0); BVoigt(4,0) = 2.0*B(2,1); BVoigt(5,0) = 2.0*B(2,0);
  BVoigt(3,0) = B(1,0); BVoigt(4,0) = B(2,1); BVoigt(5,0) = B(2,0);
  C.Multiply('N','T',ScalarAB,AVoigt,BVoigt,ScalarThis);
  C.Multiply('N','T',ScalarAB,BVoigt,AVoigt,1.0);

  // this is explicitly what the former .Multiplies do:
//  C(0,0)= ScalarThis*C(0,0) + ScalarAB * (A(0,0)*B(0,0) + B(0,0)*A(0,0));
//  C(0,1)= ScalarThis*C(0,1) + ScalarAB * (A(0,0)*B(1,1) + B(0,0)*A(1,1));
//  C(0,2)= ScalarThis*C(0,2) + ScalarAB * (A(0,0)*B(2,2) + B(0,0)*A(2,2));
//  C(0,3)= ScalarThis*C(0,3) + ScalarAB * (A(0,0)*B(1,0) + B(0,0)*A(1,0));
//  C(0,4)= ScalarThis*C(0,4) + ScalarAB * (A(0,0)*B(2,1) + B(0,0)*A(2,1));
//  C(0,5)= ScalarThis*C(0,5) + ScalarAB * (A(0,0)*B(2,0) + B(0,0)*A(2,0));
//
//  C(1,0)= ScalarThis*C(1,0) + ScalarAB * (A(1,1)*B(0,0) + B(1,1)*A(0,0));
//  C(1,1)= ScalarThis*C(1,1) + ScalarAB * (A(1,1)*B(1,1) + B(1,1)*A(1,1));
//  C(1,2)= ScalarThis*C(1,2) + ScalarAB * (A(1,1)*B(2,2) + B(1,1)*A(2,2));
//  C(1,3)= ScalarThis*C(1,3) + ScalarAB * (A(1,1)*B(1,0) + B(1,1)*A(1,0));
//  C(1,4)= ScalarThis*C(1,4) + ScalarAB * (A(1,1)*B(2,1) + B(1,1)*A(2,1));
//  C(1,5)= ScalarThis*C(1,5) + ScalarAB * (A(1,1)*B(2,0) + B(1,1)*A(2,0));
//
//  C(2,0)= ScalarThis*C(2,0) + ScalarAB * (A(2,2)*B(0,0) + B(2,2)*A(0,0));
//  C(2,1)= ScalarThis*C(2,1) + ScalarAB * (A(2,2)*B(1,1) + B(2,2)*A(1,1));
//  C(2,2)= ScalarThis*C(2,2) + ScalarAB * (A(2,2)*B(2,2) + B(2,2)*A(2,2));
//  C(2,3)= ScalarThis*C(2,3) + ScalarAB * (A(2,2)*B(1,0) + B(2,2)*A(1,0));
//  C(2,4)= ScalarThis*C(2,4) + ScalarAB * (A(2,2)*B(2,1) + B(2,2)*A(2,1));
//  C(2,5)= ScalarThis*C(2,5) + ScalarAB * (A(2,2)*B(2,0) + B(2,2)*A(2,0));
//
//  C(3,0)= ScalarThis*C(3,0) + ScalarAB * (A(1,0)*B(0,0) + B(1,0)*A(0,0));
//  C(3,1)= ScalarThis*C(3,1) + ScalarAB * (A(1,0)*B(1,1) + B(1,0)*A(1,1));
//  C(3,2)= ScalarThis*C(3,2) + ScalarAB * (A(1,0)*B(2,2) + B(1,0)*A(2,2));
//  C(3,3)= ScalarThis*C(3,3) + ScalarAB * (A(1,0)*B(1,0) + B(1,0)*A(1,0));
//  C(3,4)= ScalarThis*C(3,4) + ScalarAB * (A(1,0)*B(2,1) + B(1,0)*A(2,1));
//  C(3,5)= ScalarThis*C(3,5) + ScalarAB * (A(1,0)*B(2,0) + B(1,0)*A(2,0));
//
//  C(4,0)= ScalarThis*C(4,0) + ScalarAB * (A(2,1)*B(0,0) + B(2,1)*A(0,0));
//  C(4,1)= ScalarThis*C(4,1) + ScalarAB * (A(2,1)*B(1,1) + B(2,1)*A(1,1));
//  C(4,2)= ScalarThis*C(4,2) + ScalarAB * (A(2,1)*B(2,2) + B(2,1)*A(2,2));
//  C(4,3)= ScalarThis*C(4,3) + ScalarAB * (A(2,1)*B(1,0) + B(2,1)*A(1,0));
//  C(4,4)= ScalarThis*C(4,4) + ScalarAB * (A(2,1)*B(2,1) + B(2,1)*A(2,1));
//  C(4,5)= ScalarThis*C(4,5) + ScalarAB * (A(2,1)*B(2,0) + B(2,1)*A(2,0));
//
//  C(5,0)= ScalarThis*C(5,0) + ScalarAB * (A(2,0)*B(0,0) + B(2,0)*A(0,0));
//  C(5,1)= ScalarThis*C(5,1) + ScalarAB * (A(2,0)*B(1,1) + B(2,0)*A(1,1));
//  C(5,2)= ScalarThis*C(5,2) + ScalarAB * (A(2,0)*B(2,2) + B(2,0)*A(2,2));
//  C(5,3)= ScalarThis*C(5,3) + ScalarAB * (A(2,0)*B(1,0) + B(2,0)*A(1,0));
//  C(5,4)= ScalarThis*C(5,4) + ScalarAB * (A(2,0)*B(2,1) + B(2,0)*A(2,1));
//  C(5,5)= ScalarThis*C(5,5) + ScalarAB * (A(2,0)*B(2,0) + B(2,0)*A(2,0));

  return;
}


/*----------------------------------------------------------------------*
 | compute the "material tensor product" A o B (also known as           |
 | kronecker-tensor-product) of two 2nd order tensors                   |
 | (in matrix notation) and add the result to a 4th order tensor        |
 | (also in matrix notation) using the symmetry-conditions inherent to  |
 | material tensors, or tangent matrices, respectively                  |
 | AND the Voigt notation of E,S, and C with the famous factor 2!       |
 | The implementation is based on the Epetra-Method Matrix.Multiply.    |
 | (public)                                                    maf 11/07|
 *----------------------------------------------------------------------*/
void LINALG::ElastSymTensor_o_Multiply(Epetra_SerialDenseMatrix& C,
                                 const double ScalarAB,
                                 const Epetra_SerialDenseMatrix& A,
                                 const Epetra_SerialDenseMatrix& B,
                                 const double ScalarThis)
{
  // check sizes
  if (A.M() != A.N() || B.M() != B.N() || A.M() != 3 || B.M() != 3){
    dserror("2nd order tensors must be 3 by 3");
  }
  if (C.M() != C.N() || C.M() != 6) dserror("4th order tensor must be 6 by 6");

  /* the kronecker-product in matrix notation is:
   * A11*B11 A11*B12 A11*B13   A12*B11 A12*B12 A12*B13   A13*B11 A13*B12 A13*B13
   * A11*B21 ...
   * A11*B31 ...
   *
   * A21*B11
   * A21*B21
   * A21*B31
   * ...                                                 A33*B11 A33*B12 A33*B13
   *                                                     A33*B21 A33*B22 A33*B23
   *                                                     A33*B31 A33*B32 A33*B33
   */
  /* to reduce the resulting 9by9 matrix to 6by6 we refer to the
   * Diss. from Balzani, Anhang D, BUT
   * we consider a factor 2 for colums/rows 4-6 :
   *  C(1)               2* 1/2*(C(2)+C(3))
   *  2* 1/2*(C(2)+C(3)  2* 1/4*(C(4)+2*C(5)+C(6))
   * which is repaired later due to the "voigt-matrix":
   *    1                 1/2
   *   1/2                1/2
   */

//  C(0,0)= ScalarThis*C(0,0) + ScalarAB * A(0,0)*B(0,0);
//  C(0,1)= ScalarThis*C(0,1) + ScalarAB * A(0,0)*B(0,1);
//  C(0,2)= ScalarThis*C(0,2) + ScalarAB * A(0,0)*B(0,2);
//  C(0,3)= ScalarThis*C(0,3) + ScalarAB * (A(0,1)*B(0,0) + A(0,2)*B(0,0));
//  C(0,4)= ScalarThis*C(0,4) + ScalarAB * (A(0,1)*B(0,1) + A(0,2)*B(0,1));
//  C(0,5)= ScalarThis*C(0,5) + ScalarAB * (A(0,1)*B(0,2) + A(0,2)*B(0,2));
//
//  C(1,0)= ScalarThis*C(1,0) + ScalarAB * A(0,0)*B(1,0);
//  C(1,1)= ScalarThis*C(1,1) + ScalarAB * A(0,0)*B(1,1);
//  C(1,2)= ScalarThis*C(1,2) + ScalarAB * A(0,0)*B(1,2);
//  C(1,3)= ScalarThis*C(1,3) + ScalarAB * (A(0,1)*B(1,0) + A(0,2)*B(1,0));
//  C(1,4)= ScalarThis*C(1,4) + ScalarAB * (A(0,1)*B(1,1) + A(0,2)*B(1,1));
//  C(1,5)= ScalarThis*C(1,5) + ScalarAB * (A(0,1)*B(1,2) + A(0,2)*B(1,2));
//
//  C(2,0)= ScalarThis*C(2,0) + ScalarAB * A(0,0)*B(2,0);
//  C(2,1)= ScalarThis*C(2,1) + ScalarAB * A(0,0)*B(2,1);
//  C(2,2)= ScalarThis*C(2,2) + ScalarAB * A(0,0)*B(2,2);
//  C(2,3)= ScalarThis*C(2,3) + ScalarAB * (A(0,1)*B(2,0) + A(0,2)*B(2,0));
//  C(2,4)= ScalarThis*C(2,4) + ScalarAB * (A(0,1)*B(2,1) + A(0,2)*B(2,1));
//  C(2,5)= ScalarThis*C(2,5) + ScalarAB * (A(0,1)*B(2,2) + A(0,2)*B(2,2));
//
//  C(3,0)= ScalarThis*C(3,0) + ScalarAB * (A(1,0)*B(0,0) + A(2,0)*B(0,0));
//  C(3,1)= ScalarThis*C(3,1) + ScalarAB * (A(1,0)*B(0,1) + A(2,0)*B(0,1));
//  C(3,2)= ScalarThis*C(3,2) + ScalarAB * (A(1,0)*B(0,2) + A(2,0)*B(0,2));
//  C(3,3)= ScalarThis*C(3,3) + ScalarAB * 0.5*(A(1,1)*B(0,0) + 2.0*A(1,2)*B(0,0) + A(2,2)*B(0,0));
//  C(3,4)= ScalarThis*C(3,4) + ScalarAB * 0.5*(A(1,1)*B(0,1) + 2.0*A(1,2)*B(0,1) + A(2,2)*B(0,1));
//  C(3,5)= ScalarThis*C(3,5) + ScalarAB * 0.5*(A(1,1)*B(0,2) + 2.0*A(1,2)*B(0,2) + A(2,2)*B(0,2));
//
//  C(4,0)= ScalarThis*C(4,0) + ScalarAB * (A(1,0)*B(1,0) + A(2,0)*B(1,0));
//  C(4,1)= ScalarThis*C(4,1) + ScalarAB * (A(1,0)*B(1,1) + A(2,0)*B(1,1));
//  C(4,2)= ScalarThis*C(4,2) + ScalarAB * (A(1,0)*B(1,2) + A(2,0)*B(1,2));
//  C(4,3)= ScalarThis*C(4,3) + ScalarAB * 0.5*(A(1,1)*B(1,0) + 2.0*A(1,2)*B(1,0) + A(2,2)*B(1,0));
//  C(4,4)= ScalarThis*C(4,4) + ScalarAB * 0.5*(A(1,1)*B(1,1) + 2.0*A(1,2)*B(1,1) + A(2,2)*B(1,1));
//  C(4,5)= ScalarThis*C(4,5) + ScalarAB * 0.5*(A(1,1)*B(1,2) + 2.0*A(1,2)*B(1,2) + A(2,2)*B(1,2));
//
//  C(5,0)= ScalarThis*C(5,0) + ScalarAB * (A(1,0)*B(2,0) + A(2,0)*B(2,0));
//  C(5,1)= ScalarThis*C(5,1) + ScalarAB * (A(1,0)*B(2,1) + A(2,0)*B(2,1));
//  C(5,2)= ScalarThis*C(5,2) + ScalarAB * (A(1,0)*B(2,2) + A(2,0)*B(2,2));
//  C(5,3)= ScalarThis*C(5,3) + ScalarAB * 0.5*(A(1,1)*B(2,0) + 2.0*A(1,2)*B(2,0) + A(2,2)*B(2,0));
//  C(5,4)= ScalarThis*C(5,4) + ScalarAB * 0.5*(A(1,1)*B(2,1) + 2.0*A(1,2)*B(2,1) + A(2,2)*B(2,1));
//  C(5,5)= ScalarThis*C(5,5) + ScalarAB * 0.5*(A(1,1)*B(2,2) + 2.0*A(1,2)*B(2,2) + A(2,2)*B(2,2));

  C(0,0)= ScalarThis*C(0,0) + ScalarAB * 0.5 * (A(0,0)*B(0,0) + A(0,0)*B(0,0));//C1111
  C(0,1)= ScalarThis*C(0,1) + ScalarAB * 0.5 * (A(0,1)*B(0,1) + A(0,1)*B(0,1));//C1122
  C(0,2)= ScalarThis*C(0,2) + ScalarAB * 0.5 * (A(0,2)*B(0,2) + A(0,2)*B(0,2));//C1133
  C(0,3)= ScalarThis*C(0,3) + ScalarAB * 0.5 * (A(0,0)*B(0,1) + A(0,1)*B(0,0));//C1112
  C(0,4)= ScalarThis*C(0,4) + ScalarAB * 0.5 * (A(0,1)*B(0,2) + A(0,2)*B(0,1));//C1123
  C(0,5)= ScalarThis*C(0,5) + ScalarAB * 0.5 * (A(0,0)*B(0,2) + A(0,2)*B(0,0));//C1113

  C(1,0)= ScalarThis*C(1,0) + ScalarAB * 0.5 * (A(1,0)*B(1,0) + A(1,0)*B(1,0));//C2211
  C(1,1)= ScalarThis*C(1,1) + ScalarAB * 0.5 * (A(1,1)*B(1,1) + A(1,1)*B(1,1));//C2222
  C(1,2)= ScalarThis*C(1,2) + ScalarAB * 0.5 * (A(1,2)*B(1,2) + A(1,2)*B(1,2));//C2233
  C(1,3)= ScalarThis*C(1,3) + ScalarAB * 0.5 * (A(1,0)*B(1,1) + A(1,1)*B(1,0));//C2212
  C(1,4)= ScalarThis*C(1,4) + ScalarAB * 0.5 * (A(1,1)*B(1,2) + A(1,2)*B(1,1));//C2223
  C(1,5)= ScalarThis*C(1,5) + ScalarAB * 0.5 * (A(1,0)*B(1,2) + A(1,2)*B(1,0));//C2213

  C(2,0)= ScalarThis*C(2,0) + ScalarAB * 0.5 * (A(2,0)*B(2,0) + A(2,0)*B(2,0));//C3311
  C(2,1)= ScalarThis*C(2,1) + ScalarAB * 0.5 * (A(2,1)*B(2,1) + A(2,1)*B(2,1));//C3322
  C(2,2)= ScalarThis*C(2,2) + ScalarAB * 0.5 * (A(2,2)*B(2,2) + A(2,2)*B(2,2));//C3333
  C(2,3)= ScalarThis*C(2,3) + ScalarAB * 0.5 * (A(2,1)*B(2,1) + A(2,1)*B(2,0));//C3312
  C(2,4)= ScalarThis*C(2,4) + ScalarAB * 0.5 * (A(2,1)*B(2,2) + A(2,2)*B(2,1));//C3323
  C(2,5)= ScalarThis*C(2,5) + ScalarAB * 0.5 * (A(2,0)*B(2,2) + A(2,2)*B(2,0));//C3313

  C(3,0)= ScalarThis*C(3,0) + ScalarAB * 0.5 * (A(0,0)*B(1,0) + A(0,0)*B(1,0));//C1211
  C(3,1)= ScalarThis*C(3,1) + ScalarAB * 0.5 * (A(0,1)*B(1,1) + A(0,1)*B(1,1));//C1222
  C(3,2)= ScalarThis*C(3,2) + ScalarAB * 0.5 * (A(0,2)*B(1,2) + A(0,2)*B(1,2));//C1233
  C(3,3)= ScalarThis*C(3,3) + ScalarAB * 0.5 * (A(0,0)*B(1,1) + A(0,1)*B(1,0));//C1212
  C(3,4)= ScalarThis*C(3,4) + ScalarAB * 0.5 * (A(0,1)*B(1,2) + A(0,2)*B(1,1));//C1223
  C(3,5)= ScalarThis*C(3,5) + ScalarAB * 0.5 * (A(0,0)*B(1,2) + A(0,2)*B(1,0));//C1213

  C(4,0)= ScalarThis*C(4,0) + ScalarAB * 0.5 * (A(1,0)*B(2,0) + A(1,0)*B(2,0));//C2311
  C(4,1)= ScalarThis*C(4,1) + ScalarAB * 0.5 * (A(1,1)*B(2,1) + A(1,1)*B(2,1));//C2322
  C(4,2)= ScalarThis*C(4,2) + ScalarAB * 0.5 * (A(1,2)*B(2,2) + A(1,2)*B(2,2));//C2333
  C(4,3)= ScalarThis*C(4,3) + ScalarAB * 0.5 * (A(1,0)*B(2,1) + A(1,1)*B(2,0));//C2312
  C(4,4)= ScalarThis*C(4,4) + ScalarAB * 0.5 * (A(1,1)*B(2,2) + A(1,2)*B(2,1));//C2323
  C(4,5)= ScalarThis*C(4,5) + ScalarAB * 0.5 * (A(1,0)*B(2,2) + A(1,2)*B(2,0));//C2313

  C(5,0)= ScalarThis*C(5,0) + ScalarAB * 0.5 * (A(0,0)*B(2,0) + A(0,0)*B(2,0));//C1311
  C(5,1)= ScalarThis*C(5,1) + ScalarAB * 0.5 * (A(0,1)*B(2,1) + A(0,1)*B(2,1));//C1322
  C(5,2)= ScalarThis*C(5,2) + ScalarAB * 0.5 * (A(0,2)*B(2,2) + A(0,2)*B(2,2));//C1333
  C(5,3)= ScalarThis*C(5,3) + ScalarAB * 0.5 * (A(0,0)*B(2,1) + A(0,1)*B(2,0));//C1312
  C(5,4)= ScalarThis*C(5,4) + ScalarAB * 0.5 * (A(0,1)*B(2,2) + A(0,2)*B(2,1));//C1323
  C(5,5)= ScalarThis*C(5,5) + ScalarAB * 0.5 * (A(0,0)*B(2,2) + A(0,2)*B(2,0));//C1313

  return;

}


/*----------------------------------------------------------------------*
| invert a dense nonsymmetric matrix (public)       g.bau 03/07|
*----------------------------------------------------------------------*/
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
void LINALG::ApplyDirichlettoSystem(RefCountPtr<Epetra_Vector>&      x,
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
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(RCP<LINALG::SparseOperator>      A,
                                    RefCountPtr<Epetra_Vector>&      x,
                                    RefCountPtr<Epetra_Vector>&      b,
                                    const RefCountPtr<Epetra_Vector> dbcval,
                                    const RefCountPtr<Epetra_Vector> dbctoggle)
{
  A->ApplyDirichlet(dbctoggle);
  ApplyDirichlettoSystem(x,b,dbcval,dbctoggle);
}


/*----------------------------------------------------------------------*
 | split matrix into 2x2 block system                              06/06|
 *----------------------------------------------------------------------*/
bool LINALG::SplitMatrix2x2(RCP<Epetra_CrsMatrix> A,
                            RCP<Epetra_Map>& A11rowmap,
                            RCP<Epetra_Map>& A22rowmap,
                            RCP<Epetra_CrsMatrix>& A11,
                            RCP<Epetra_CrsMatrix>& A12,
                            RCP<Epetra_CrsMatrix>& A21,
                            RCP<Epetra_CrsMatrix>& A22)
{
  if (A==null)
    dserror("LINALG::SplitMatrix2x2: A==null on entry");

  if (A11rowmap==null && A22rowmap != null)
    A11rowmap = LINALG::SplitMap(A->RowMap(),*A22rowmap);
  else if (A11rowmap != null && A22rowmap != null);
  else if (A11rowmap != null && A22rowmap == null)
    A22rowmap = LINALG::SplitMap(A->RowMap(),*A11rowmap);
  else
  	dserror("LINALG::SplitMatrix2x2: Both A11rowmap and A22rowmap == null on entry");

  const Epetra_Comm& Comm   = A->Comm();
  const Epetra_Map&  A22map = *(A22rowmap.get());
  const Epetra_Map&  A11map = *(A11rowmap.get());

  //----------------------------- create a parallel redundant map of A22map
  map<int,int> a22gmap;
  {
    vector<int> a22global(A22map.NumGlobalElements());
    int count=0;
    for (int proc=0; proc<Comm.NumProc(); ++proc)
    {
      int length = 0;
      if (proc==Comm.MyPID())
      {
        for (int i=0; i<A22map.NumMyElements(); ++i)
        {
          a22global[count+length] = A22map.GID(i);
          ++length;
        }
      }
      Comm.Broadcast(&length,1,proc);
      Comm.Broadcast(&a22global[count],length,proc);
      count += length;
    }
    if (count != A22map.NumGlobalElements())
    	dserror("LINALG::SplitMatrix2x2: mismatch in dimensions");

    // create the map
    for (int i=0; i<count; ++i)
      a22gmap[a22global[i]] = 1;
    a22global.clear();
  }

  //--------------------------------------------------- create matrix A22
  A22 = rcp(new Epetra_CrsMatrix(Copy,A22map,100));
  {
    vector<int>    a22gcindices(100);
    vector<double> a22values(100);
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A22map.MyGID(grid)==false)
        continue;
      //cout << "Row " << grid << " in A22 Columns ";
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->ExtractMyRowView(i,numentries,values,cindices);
      if (err)
      	dserror("LINALG::SplitMatrix2x2: A->ExtractMyRowView returned &i",err);

      if (numentries>(int)a22gcindices.size())
      {
        a22gcindices.resize(numentries);
        a22values.resize(numentries);
      }
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid in a22gmap
        map<int,int>::iterator curr = a22gmap.find(gcid);
        if (curr==a22gmap.end()) continue;
        //cout << gcid << " ";
        a22gcindices[count] = gcid;
        a22values[count]    = values[j];
        ++count;
      }
      //cout << endl; fflush(stdout);
      // add this filtered row to A22
      err = A22->InsertGlobalValues(grid,count,&a22values[0],&a22gcindices[0]);
      if (err<0)
      	dserror("LINALG::SplitMatrix2x2: A->InsertGlobalValues returned &i",err);

    } //for (int i=0; i<A->NumMyRows(); ++i)
    a22gcindices.clear();
    a22values.clear();
  }
  A22->FillComplete();
  A22->OptimizeStorage();

  //----------------------------------------------------- create matrix A11
  A11 = rcp(new Epetra_CrsMatrix(Copy,A11map,100));
  {
    vector<int>    a11gcindices(100);
    vector<double> a11values(100);
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A11map.MyGID(grid)==false) continue;
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->ExtractMyRowView(i,numentries,values,cindices);
      if (err)
      	dserror("LINALG::SplitMatrix2x2: A->ExtractMyRowView returned &i",err);

      if (numentries>(int)a11gcindices.size())
      {
        a11gcindices.resize(numentries);
        a11values.resize(numentries);
      }
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid as part of a22gmap
        map<int,int>::iterator curr = a22gmap.find(gcid);
        if (curr!=a22gmap.end()) continue;
        a11gcindices[count] = gcid;
        a11values[count] = values[j];
        ++count;
      }
      err = A11->InsertGlobalValues(grid,count,&a11values[0],&a11gcindices[0]);
      if (err<0)
      	dserror("LINALG::SplitMatrix2x2: A->InsertGlobalValues returned &i",err);

    } // for (int i=0; i<A->NumMyRows(); ++i)
    a11gcindices.clear();
    a11values.clear();
  }
  A11->FillComplete();
  A11->OptimizeStorage();

  //---------------------------------------------------- create matrix A12
  A12 = rcp(new Epetra_CrsMatrix(Copy,A11map,100));
  {
    vector<int>    a12gcindices(100);
    vector<double> a12values(100);
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A11map.MyGID(grid)==false) continue;
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->ExtractMyRowView(i,numentries,values,cindices);
      if (err)
      	dserror("LINALG::SplitMatrix2x2: A->ExtractMyRowView returned &i",err);

      if (numentries>(int)a12gcindices.size())
      {
        a12gcindices.resize(numentries);
        a12values.resize(numentries);
      }
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid as part of a22gmap
        map<int,int>::iterator curr = a22gmap.find(gcid);
        if (curr==a22gmap.end()) continue;
        a12gcindices[count] = gcid;
        a12values[count] = values[j];
        ++count;
      }
      err = A12->InsertGlobalValues(grid,count,&a12values[0],&a12gcindices[0]);
      if (err<0)
      	dserror("LINALG::SplitMatrix2x2: A->InsertGlobalValues returned &i",err);

    } // for (int i=0; i<A->NumMyRows(); ++i)
    a12values.clear();
    a12gcindices.clear();
  }
  A12->FillComplete(A22map,A11map);
  A12->OptimizeStorage();

  //----------------------------------------------------------- create A21
  A21 = rcp(new Epetra_CrsMatrix(Copy,A22map,100));
  {
    vector<int>    a21gcindices(100);
    vector<double> a21values(100);
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A22map.MyGID(grid)==false) continue;
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->ExtractMyRowView(i,numentries,values,cindices);
      if (err)
      	dserror("LINALG::SplitMatrix2x2: A->ExtractMyRowView returned &i",err);

      if (numentries>(int)a21gcindices.size())
      {
        a21gcindices.resize(numentries);
        a21values.resize(numentries);
      }
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid as part of a22gmap
        map<int,int>::iterator curr = a22gmap.find(gcid);
        if (curr!=a22gmap.end()) continue;
        a21gcindices[count] = gcid;
        a21values[count] = values[j];
        ++count;
      }
      err = A21->InsertGlobalValues(grid,count,&a21values[0],&a21gcindices[0]);
      if (err<0)
      	dserror("LINALG::SplitMatrix2x2: A->InsertGlobalValues returned &i",err);

    } // for (int i=0; i<A->NumMyRows(); ++i)
    a21values.clear();
    a21gcindices.clear();
  }
  A21->FillComplete(A11map,A22map);
  A21->OptimizeStorage();

  //-------------------------------------------------------------- tidy up
  a22gmap.clear();
  return true;
}

/*----------------------------------------------------------------------*
 | split matrix into 2x2 block system                         popp 02/08|
 *----------------------------------------------------------------------*/
bool LINALG::SplitMatrix2x2(RCP<Epetra_CrsMatrix> A,
                            RCP<Epetra_Map>& A11rowmap,
                            RCP<Epetra_Map>& A22rowmap,
                            RCP<Epetra_Map>& A11domainmap,
                            RCP<Epetra_Map>& A22domainmap,
                            RCP<Epetra_CrsMatrix>& A11,
                            RCP<Epetra_CrsMatrix>& A12,
                            RCP<Epetra_CrsMatrix>& A21,
                            RCP<Epetra_CrsMatrix>& A22)
{
  if (A==null)
    dserror("LINALG::SplitMatrix2x2: A==null on entry");

  // check and complete input row maps
  if (A11rowmap==null && A22rowmap != null)
    A11rowmap = LINALG::SplitMap(A->RowMap(),*A22rowmap);
  else if (A11rowmap != null && A22rowmap != null);
  else if (A11rowmap != null && A22rowmap == null)
    A22rowmap = LINALG::SplitMap(A->RowMap(),*A11rowmap);
  else
  	dserror("LINALG::SplitMatrix2x2: Both A11rowmap and A22rowmap == null on entry");

  // check and complete input domain maps
  if (A11domainmap==null && A22domainmap != null)
  	A11domainmap = LINALG::SplitMap(A->DomainMap(),*A22domainmap);
  else if (A11domainmap != null && A22domainmap != null);
  else if (A11domainmap != null && A22domainmap == null)
    A22domainmap = LINALG::SplitMap(A->DomainMap(),*A11domainmap);
  else
  	dserror("LINALG::SplitMatrix2x2: Both A11domainmap and A22domainmap == null on entry");

  // local variables
  const Epetra_Comm& Comm   = A->Comm();
  const Epetra_Map&  A11rmap = *(A11rowmap.get());
  const Epetra_Map&  A11dmap = *(A11domainmap.get());
  const Epetra_Map&  A22rmap = *(A22rowmap.get());
  const Epetra_Map&  A22dmap = *(A22domainmap.get());

  //----------------------------- create a parallel redundant map of A11domainmap
  map<int,int> a11gmap;
  {
    vector<int> a11global(A11dmap.NumGlobalElements());
    int count=0;
    for (int proc=0; proc<Comm.NumProc(); ++proc)
    {
      int length = 0;
      if (proc==Comm.MyPID())
      {
        for (int i=0; i<A11dmap.NumMyElements(); ++i)
        {
          a11global[count+length] = A11dmap.GID(i);
          ++length;
        }
      }
      Comm.Broadcast(&length,1,proc);
      Comm.Broadcast(&a11global[count],length,proc);
      count += length;
    }
    if (count != A11dmap.NumGlobalElements())
    	dserror("LINALG::SplitMatrix2x2: mismatch in dimensions");

    // create the map
    for (int i=0; i<count; ++i)
      a11gmap[a11global[i]] = 1;
    a11global.clear();
  }

  //----------------------------------------------------- create matrix A11
  if (A11rmap.NumGlobalElements()>0 && A11dmap.NumGlobalElements()>0)
  {
  	A11 = rcp(new Epetra_CrsMatrix(Copy,A11rmap,100));
  	{
  		vector<int>    a11gcindices(100);
  		vector<double> a11values(100);
  		for (int i=0; i<A->NumMyRows(); ++i)
  		{
  			const int grid = A->GRID(i);
  			if (A11rmap.MyGID(grid)==false) continue;
  			int     numentries;
  			double* values;
  			int*    cindices;
  			int err = A->ExtractMyRowView(i,numentries,values,cindices);
  			if (err)
  				dserror("LINALG::SplitMatrix2x2: A->ExtractMyRowView returned %i",err);

  			if (numentries>(int)a11gcindices.size())
  			{
  				a11gcindices.resize(numentries);
  				a11values.resize(numentries);
  			}
  			int count=0;
  			for (int j=0; j<numentries; ++j)
  			{
  				const int gcid = A->ColMap().GID(cindices[j]);
  				// see whether we have gcid as part of a11gmap
  				map<int,int>::iterator curr = a11gmap.find(gcid);
  				if (curr==a11gmap.end()) continue;
  				a11gcindices[count] = gcid;
  				a11values[count] = values[j];
  				++count;
  			}
  			err = A11->InsertGlobalValues(grid,count,&a11values[0],&a11gcindices[0]);
  			if (err<0)
  				dserror("LINALG::SplitMatrix2x2: A->InsertGlobalValues returned %i",err);

  		} // for (int i=0; i<A->NumMyRows(); ++i)
  		a11gcindices.clear();
  		a11values.clear();
  	}
  	LINALG::Complete(*A11,A11dmap,A11rmap);
  }

  //--------------------------------------------------- create matrix A22
  if (A22rmap.NumGlobalElements()>0 && A22dmap.NumGlobalElements()>0)
  {
	  A22 = rcp(new Epetra_CrsMatrix(Copy,A22rmap,100));
	  {
	    vector<int>    a22gcindices(100);
	    vector<double> a22values(100);
	    for (int i=0; i<A->NumMyRows(); ++i)
	    {
	      const int grid = A->GRID(i);
	      if (A22rmap.MyGID(grid)==false) continue;
	      int     numentries;
	      double* values;
	      int*    cindices;
	      int err = A->ExtractMyRowView(i,numentries,values,cindices);
	      if (err)
	      	dserror("LINALG::SplitMatrix2x2: A->ExtractMyRowView returned %i",err);

	      if (numentries>(int)a22gcindices.size())
	      {
	        a22gcindices.resize(numentries);
	        a22values.resize(numentries);
	      }
	      int count=0;
	      for (int j=0; j<numentries; ++j)
	      {
	        const int gcid = A->ColMap().GID(cindices[j]);
	        // see whether we have gcid as part of a11gmap
	        map<int,int>::iterator curr = a11gmap.find(gcid);
	        if (curr!=a11gmap.end()) continue;
	        a22gcindices[count] = gcid;
	        a22values[count]    = values[j];
	        ++count;
	      }
	      err = A22->InsertGlobalValues(grid,count,&a22values[0],&a22gcindices[0]);
	      if (err<0)
	      	dserror("LINALG::SplitMatrix2x2: A->InsertGlobalValues returned %i",err);

	    } //for (int i=0; i<A->NumMyRows(); ++i)
	    a22gcindices.clear();
	    a22values.clear();
	  }
	  LINALG::Complete(*A22,A22dmap,A22rmap);
  }

  //---------------------------------------------------- create matrix A12
  if (A11rmap.NumGlobalElements()>0 && A22dmap.NumGlobalElements()>0)
  {
	  A12 = rcp(new Epetra_CrsMatrix(Copy,A11rmap,100));
	  {
	    vector<int>    a12gcindices(100);
	    vector<double> a12values(100);
	    for (int i=0; i<A->NumMyRows(); ++i)
	    {
	      const int grid = A->GRID(i);
	      if (A11rmap.MyGID(grid)==false) continue;
	      int     numentries;
	      double* values;
	      int*    cindices;
	      int err = A->ExtractMyRowView(i,numentries,values,cindices);
	      if (err)
	      	dserror("LINALG::SplitMatrix2x2: A->ExtractMyRowView returned %i",err);

	      if (numentries>(int)a12gcindices.size())
	      {
	        a12gcindices.resize(numentries);
	        a12values.resize(numentries);
	      }
	      int count=0;
	      for (int j=0; j<numentries; ++j)
	      {
	        const int gcid = A->ColMap().GID(cindices[j]);
	        // see whether we have gcid as part of a11gmap
	        map<int,int>::iterator curr = a11gmap.find(gcid);
	        if (curr!=a11gmap.end()) continue;
	        a12gcindices[count] = gcid;
	        a12values[count] = values[j];
	        ++count;
	      }
	      err = A12->InsertGlobalValues(grid,count,&a12values[0],&a12gcindices[0]);
	      if (err<0)
	      	dserror("LINALG::SplitMatrix2x2: A->InsertGlobalValues returned %i",err);

	    } // for (int i=0; i<A->NumMyRows(); ++i)
	    a12values.clear();
	    a12gcindices.clear();
	  }
	  LINALG::Complete(*A12,A22dmap,A11rmap);
  }

  //----------------------------------------------------------- create A21
  if (A22rmap.NumGlobalElements()>0 && A11dmap.NumGlobalElements()>0)
  {
	  A21 = rcp(new Epetra_CrsMatrix(Copy,A22rmap,100));
	  {
	    vector<int>    a21gcindices(100);
	    vector<double> a21values(100);
	    for (int i=0; i<A->NumMyRows(); ++i)
	    {
	      const int grid = A->GRID(i);
	      if (A22rmap.MyGID(grid)==false) continue;
	      int     numentries;
	      double* values;
	      int*    cindices;
	      int err = A->ExtractMyRowView(i,numentries,values,cindices);
	      if (err)
	      	dserror("LINALG::SplitMatrix2x2: A->ExtractMyRowView returned %i",err);

	      if (numentries>(int)a21gcindices.size())
	      {
	        a21gcindices.resize(numentries);
	        a21values.resize(numentries);
	      }
	      int count=0;
	      for (int j=0; j<numentries; ++j)
	      {
	        const int gcid = A->ColMap().GID(cindices[j]);
	        // see whether we have gcid as part of a11gmap
	        map<int,int>::iterator curr = a11gmap.find(gcid);
	        if (curr==a11gmap.end()) continue;
	        a21gcindices[count] = gcid;
	        a21values[count] = values[j];
	        ++count;
	      }
	      err = A21->InsertGlobalValues(grid,count,&a21values[0],&a21gcindices[0]);
	      if (err<0)
	      	dserror("LINALG::SplitMatrix2x2: A->InsertGlobalValues returned %i",err);

	    } // for (int i=0; i<A->NumMyRows(); ++i)
	    a21values.clear();
	    a21gcindices.clear();
	  }
	  LINALG::Complete(*A21,A11dmap,A22rmap);
  }

  //-------------------------------------------------------------- tidy up
  a11gmap.clear();
  return true;
}

/*----------------------------------------------------------------------*
 | split a map into 2 pieces with given Agiven                     06/06|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::SplitMap(const Epetra_Map& Amap,
                                          const Epetra_Map& Agiven)
{
  const Epetra_Comm& Comm = Amap.Comm();
  const Epetra_Map&  Ag = Agiven;

  int count=0;
  vector<int> myaugids(Amap.NumMyElements());
  for (int i=0; i<Amap.NumMyElements(); ++i)
  {
    const int gid = Amap.GID(i);
    if (Ag.MyGID(gid)) continue;
    myaugids[count] = gid;
    ++count;
  }
  myaugids.resize(count);
  int gcount;
  Comm.SumAll(&count,&gcount,1);
  Teuchos::RCP<Epetra_Map> Aunknown = Teuchos::rcp(new Epetra_Map(gcount,count,&myaugids[0],0,Comm));
  myaugids.clear();
  return Aunknown;
}


/*----------------------------------------------------------------------*
 | merge two given maps to one map                            popp 01/08|
 *----------------------------------------------------------------------*/
RCP<Epetra_Map> LINALG::MergeMap(const Epetra_Map& map1,
                                 const Epetra_Map& map2)
{
  // check for unique GIDs and for identity
  if ((!map1.UniqueGIDs()) || (!map2.UniqueGIDs()))
    dserror("LINALG::MergeMap: One or both input maps are not unique");
  if (map1.SameAs(map2))
    return rcp(new Epetra_Map(map1));

  vector<int> mygids(map1.NumMyElements()+map2.NumMyElements());
  int count = map1.NumMyElements();

  // get GIDs of input map1
  for (int i=0;i<count;++i)
    mygids[i] = map1.GID(i);

  // add GIDs of input map2 (only new ones)
  for (int i=0;i<map2.NumMyElements();++i)
    if (!map1.MyGID(map2.GID(i)))
    {
      mygids[count]=map2.GID(i);
      ++count;
    }
  mygids.resize(count);

	// sort merged map
	sort(mygids.begin(),mygids.end());

	return rcp(new Epetra_Map(-1,(int)mygids.size(),&mygids[0],0,map1.Comm()));
}

/*----------------------------------------------------------------------*
 | merge two given maps to one map                            popp 01/08|
 *----------------------------------------------------------------------*/
RCP<Epetra_Map> LINALG::MergeMap(const RCP<Epetra_Map>& map1,
                                 const RCP<Epetra_Map>& map2)
{
  // check for cases with null RCPs
  if (map1==null && map2==null)
    return null;
  else if (map1==null)
    return rcp(new Epetra_Map(*map2));
  else if (map2==null)
    return rcp(new Epetra_Map(*map1));

  // wrapped call to non-RCP version of MergeMap
  return LINALG::MergeMap(*map1,*map2);
}

/*----------------------------------------------------------------------*
 | split a vector into 2 pieces with given submaps                 06/06|
 *----------------------------------------------------------------------*/
bool LINALG::SplitVector(const Epetra_Vector& x,
                         const Epetra_Map& x1map,
                         Epetra_Vector*&   x1,
                         const Epetra_Map& x2map,
                         Epetra_Vector*&   x2)
{
  x1 = new Epetra_Vector(x1map,false);
  x2 = new Epetra_Vector(x2map,false);

  //use an exporter or importer object
  Epetra_Export exporter_x1(x.Map(),x1map);
  Epetra_Export exporter_x2(x.Map(),x2map);

  int err = x1->Export(x,exporter_x1,Insert);
  if (err)
  {
    cout << "***ERR*** LINALG::SplitVector:\n"
         << "***ERR*** Export returned " << err << endl
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  err = x2->Export(x,exporter_x2,Insert);
  if (err)
  {
    cout << "***ERR*** LINALG::SplitVector:\n"
         << "***ERR*** Export returned " << err << endl
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::PrintSparsityToPostscript(const Epetra_RowMatrix& A)
{
  Ifpack_PrintSparsity(A);
  return;
}




/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int LINALG::FindMyPos(int nummyelements, const Epetra_Comm& comm)
{
  const int myrank  = comm.MyPID();
  const int numproc = comm.NumProc();

  vector<int> snum(numproc);
  vector<int> rnum(numproc);
  fill(snum.begin(), snum.end(), 0);
  snum[myrank] = nummyelements;

  comm.SumAll(&snum[0],&rnum[0],numproc);

  return std::accumulate(&rnum[0], &rnum[myrank], 0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::AllreduceEMap(vector<int>& rredundant, const Epetra_Map& emap)
{
  int mynodepos = FindMyPos(emap.NumMyElements(), emap.Comm());

  vector<int> sredundant(emap.NumGlobalElements());
  fill(sredundant.begin(), sredundant.end(), 0);

  int* gids = emap.MyGlobalElements();
  copy(gids, gids+emap.NumMyElements(), &sredundant[mynodepos]);

  rredundant.resize(emap.NumGlobalElements());
  emap.Comm().SumAll(&sredundant[0], &rredundant[0], emap.NumGlobalElements());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::AllreduceEMap(map<int,int>& idxmap, const Epetra_Map& emap)
{
  idxmap.clear();

  vector<int> rredundant(emap.NumGlobalElements());
  AllreduceEMap(rredundant, emap);

  for (unsigned i=0; i<rredundant.size(); ++i)
  {
    idxmap[rredundant[i]] = i;
  }
}

/*----------------------------------------------------------------------*
 |  create an allreduced map on a distinct processor (public)  gjb 12/07|
 *----------------------------------------------------------------------*/
RCP<Epetra_Map> LINALG::AllreduceEMap(const Epetra_Map& emap, const int pid)
{
  vector<int> rv;
  AllreduceEMap(rv,emap);
  RefCountPtr<Epetra_Map> rmap;

  if (emap.Comm().MyPID()==pid)
  {
	  rmap = rcp(new Epetra_Map(-1,rv.size(),&rv[0],0,emap.Comm()));
	  // check the map
	  dsassert(rmap->NumMyElements() == rmap->NumGlobalElements(),
	  			  "Processor with pid does not get all map elements");
  }
  else
  {
	  rv.clear();
	  rmap = rcp(new Epetra_Map(-1,0,NULL,0,emap.Comm()));
	  // check the map
	  dsassert(rmap->NumMyElements() == 0,
	  			  "At least one proc will keep a map element");
  }
  return rmap;
}


/*----------------------------------------------------------------------*
 |  Send and receive lists of ints.  (heiner 09/07)                     |
 *----------------------------------------------------------------------*/
void LINALG::AllToAllCommunication( const Epetra_Comm& comm,
                                    const vector< vector<int> >& send,
                                    vector< vector<int> >& recv )
{
#ifndef PARALLEL

  dsassert(send.size()==1, "there has to be just one entry for sending");

  // make a copy
  recv.clear();
  recv.push_back(send[0]);

#else

  if (comm.NumProc()==1)
  {
    dsassert(send.size()==1, "there has to be just one entry for sending");

    // make a copy
    recv.clear();
    recv.push_back(send[0]);
  }
  else
  {
    const Epetra_MpiComm& mpicomm = dynamic_cast<const Epetra_MpiComm&>(comm);

    vector<int> sendbuf;
    vector<int> sendcounts;
    sendcounts.reserve( comm.NumProc() );
    vector<int> sdispls;
    sdispls.reserve( comm.NumProc() );

    int displacement = 0;
    sdispls.push_back( 0 );
    for ( vector< vector<int> >::const_iterator iter = send.begin();
          iter != send.end(); ++iter )
    {
        sendbuf.insert( sendbuf.end(), iter->begin(), iter->end() );
        sendcounts.push_back( iter->size() );
        displacement += iter->size();
        sdispls.push_back( displacement );
    }

    vector<int> recvcounts( comm.NumProc() );

    // initial communication: Request. Send and receive the number of
    // ints we communicate with each process.

    int status = MPI_Alltoall( &sendcounts[0], 1, MPI_INT,
                               &recvcounts[0], 1, MPI_INT, mpicomm.GetMpiComm() );

    if ( status != MPI_SUCCESS )
        dserror( "MPI_Alltoall returned status=%d", status );

    vector<int> rdispls;
    rdispls.reserve( comm.NumProc() );

    displacement = 0;
    rdispls.push_back( 0 );
    for ( vector<int>::const_iterator iter = recvcounts.begin();
          iter != recvcounts.end(); ++iter )
    {
        displacement += *iter;
        rdispls.push_back( displacement );
    }

    vector<int> recvbuf( rdispls.back() );

    // transmit communication: Send and get the data.

    status = MPI_Alltoallv ( &sendbuf[0], &sendcounts[0], &sdispls[0], MPI_INT,
                             &recvbuf[0], &recvcounts[0], &rdispls[0], MPI_INT,
                             mpicomm.GetMpiComm() );
    if ( status != MPI_SUCCESS )
        dserror( "MPI_Alltoallv returned status=%d", status );

    recv.clear();
    for ( int proc = 0; proc < comm.NumProc(); ++proc )
    {
        recv.push_back( vector<int>( &recvbuf[rdispls[proc]], &recvbuf[rdispls[proc+1]] ) );
    }
  }

#endif // PARALLEL
}

#endif  // #ifdef CCADISCRET
