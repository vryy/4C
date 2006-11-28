/*!----------------------------------------------------------------------
\file design.cpp
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

#include "design.H"
#include "designdiscretization.H"
#include "dserror.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Design::Design(RefCountPtr<Epetra_Comm> comm) :
comm_(comm)
{
  for (int i=0; i<3; ++i) entity_[i] = rcp(new DesignDiscretization(comm_));
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Design::Design(const DRT::Design& old)
{
  dserror("Do not have copy-ctor");
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Design::~Design()
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::DesignDiscretization> DRT::Design::operator [] (const int index) const
{ 
  if (index != 0 && index != 1 && index != 2)
    dserror("index out of range, has to be 0 for lines/nodes, 1 for surfaces, 2 for volumes"); 
  return entity_[index]; 
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
