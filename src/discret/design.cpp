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
#include "dserror.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Design::Design(RefCountPtr<Epetra_Comm> comm) :
comm_(comm)
{
  for (int i=0; i<3; ++i) entity_[i] = rcp(new DesignDiscretization(comm_));
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Design::Design(const CCADISCRETIZATION::Design& old) :
comm_(old.comm_)
{
  // make a deep copy, do not share data between old and this
  for (int i=0; i<3; ++i)
    entity_[i] = rcp(new DesignDiscretization(*(old.entity_[i])));
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Design::~Design()
{
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
