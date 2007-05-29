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
      //if (abs(val)<1.0e-10) continue; // do not assemble zeros
      if (abs(val)==0) continue; // do not assemble zeros
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




/*----------------------------------------------------------------------*
| invert a dense nonsymmetric matrix (public)       g.bau 03/07|
*----------------------------------------------------------------------*/
#include <Epetra_SerialDenseSolver.h>

void LINALG::NonSymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim)
{
  if (A.M() != A.N()) dserror("Matrix is not square");
  if (A.M() != dim) dserror("Dimension supplied does not match matrix");


  //------------------------------------------------------------------------------
// this routine has to be reviewed and beautified !!!!!!!!!!!!!!!
// in the current version exactly the same results as in the old
// discretization are obtained.
// => necessary for debugging of fluid3 element.
//------------------------------------------------------------------------------
//  Take care of:   A^{-1} =  A^T^{-1}^T
//  But in numerics there will be differences, of course.



 // more elegant solution to the ordering problem mentioned below
 // could be:
/*
 Epetra_SerialDenseSolver solver;
 solver.SetMatrix(A);
int err = solver.Invert();
 if (err!=0)
  dserror("Inversion of nonsymmetric matrix failed.");
*/

  // within Epetra_SerialDenseMatrix entries are stored columnwise. Therefore
  // it is necessary to transpose before and after the LAPACK routines
  // to have the same order (rowwise!!) as in the old C code (see src/math/math1.c)

  #if 1
  double mywork=0.0;
  for (int i=0; i<dim; ++i)
    {
    for (int j=0; j<i; ++j)
    {mywork = A(i,j);
      A(i,j)=A(j,i);
      A(j,i)=mywork;
     }
    }
   #endif

  double* a = A.A();
  char uplo[5]; strncpy(uplo,"L ",2);
  vector<int> ipiv(dim);
  int lwork = 6;
  //int lwork = 10*dim;
  vector<double> work(lwork);
  int info=0;
  int n = dim;
  int m = dim;


  #if 0
  for (int i=0;i<36;++i)
  printf("a[%d] %22.16e\n",i,a[i]);
   #endif

  dgetrf(&m,&n,a,&m,&(ipiv[0]),&info);
  if (info) dserror("dgetrf returned info=%d",info);

  dgetri(&n,a,&n,&(ipiv[0]),&(work[0]),&lwork,&info);
  if (info) dserror("dgetri returned info=%d",info);

  // undo transpose of A
   #if 1
  mywork=0.0;
  for (int i=0; i<dim; ++i)
    {
    for (int j=0; j<i; ++j)
    {mywork = A(i,j);
      A(i,j)=A(j,i);
      A(j,i)=mywork;
     }
    }
   #endif

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
void LINALG::ApplyDirichlettoSystem(RefCountPtr<Epetra_CrsMatrix>&   A,
                                    RefCountPtr<Epetra_Vector>&      x,
                                    RefCountPtr<Epetra_Vector>&      b,
                                    const RefCountPtr<Epetra_Vector> dbcval,
                                    const RefCountPtr<Epetra_Vector> dbctoggle)
{
  const Epetra_Map& rowmap = A->RowMap();
#ifdef DEBUG
  if (x != null)
    if (!rowmap.PointSameAs(x->Map())) dserror("x does not match A");
  if (b != null)
    if (!rowmap.PointSameAs(b->Map())) dserror("b does not match A");
  if (dbcval != null)
    if (!rowmap.PointSameAs(dbcval->Map())) dserror("dbcval does not match A");
  if (!rowmap.PointSameAs(dbctoggle->Map())) dserror("dbctoggle does not match A");
  if (A->Filled()!=true) dserror("FillComplete was not called on A");
#endif

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

  // allocate a new matrix and copy all rows that are not dirichlet
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

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
