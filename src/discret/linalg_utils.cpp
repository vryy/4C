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
