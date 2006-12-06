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



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
