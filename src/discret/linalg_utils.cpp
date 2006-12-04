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



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
