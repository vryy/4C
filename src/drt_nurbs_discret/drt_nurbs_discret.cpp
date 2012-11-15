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


#include "drt_nurbs_discret.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::NurbsDiscretization::NurbsDiscretization(
  const string             name,
  RCP<Epetra_Comm> comm)
  :
  DRT::Discretization::Discretization(name,comm    ),
  npatches_                          (            0),
  knots_                             (Teuchos::null)
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
(RCP<DRT::NURBS::Knotvector> knots)
{

  if(knots==Teuchos::null)
  {
    dserror("trying to set invalid knotvector (%s)\n",(this->Name()).c_str());
  }

  knots_=knots;
  return;
}

/*----------------------------------------------------------------------*
 |  get a pointer to knotvector from the discretization         (public)|
 |                                                           gammi 05/08|
 *----------------------------------------------------------------------*/
RCP<DRT::NURBS::Knotvector>
DRT::NURBS::NurbsDiscretization::GetKnotVector
()
{
  if(knots_==Teuchos::null)
  {
    dserror("knotvector invalid (%s)\n",(this->Name()).c_str());
  }
  return knots_;
}

