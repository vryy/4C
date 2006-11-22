/*!----------------------------------------------------------------------
\file discret_metis.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "discret.H"
#include "exporter.H"
#include "dserror.H"




/*----------------------------------------------------------------------*
 |  redistribute discretization using metis (public)         mwgee 11/06|
 *----------------------------------------------------------------------*/
int CCADISCRETIZATION::Discretization::DistributeUsingMetis(
                                   RefCountPtr<Epetra_Map>& newnoderowmap,
                                   RefCountPtr<Epetra_Map>& newelerowmap)
{
  if (!Filled()) dserror("FillComplete() was not called on this discretization");

  // we get everything on proc 0 here to do serial metis  
  
  
  
  return 0;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
