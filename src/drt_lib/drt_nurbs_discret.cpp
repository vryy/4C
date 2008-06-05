/*!----------------------------------------------------------------------
\file drt_nurbs_discret.cpp

\brief a class to manage one nurbs discretization

<pre>
Maintainer: Peter Gamnitzer
            gammi@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_nurbs_discret.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::NurbsDiscretization::NurbsDiscretization(
  const string             name, 
  RefCountPtr<Epetra_Comm> comm) 
  :
  DRT::Discretization::Discretization(name,comm)
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::NurbsDiscretization::~NurbsDiscretization()
{
  return;
}


/*----------------------------------------------------------------------*
 |  add a knotvector to the discretization (public)          gammi 05/08|
 *----------------------------------------------------------------------*/
void 
DRT::NURBS::NurbsDiscretization::SetKnotVector
(RefCountPtr<DRT::NURBS::Knotvector> knots)
{
  knots_=knots;
  return;
}

/*----------------------------------------------------------------------*
 |  get a pointer to knotvector from the discretization         (public)|
 |                                                           gammi 05/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::NURBS::Knotvector> 
DRT::NURBS::NurbsDiscretization::GetKnotVector
()
{
  return knots_;
}

#endif
