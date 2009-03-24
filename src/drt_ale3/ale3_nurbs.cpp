/*!----------------------------------------------------------------------
\file ale3_nurbs.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#ifdef CCADISCRET

#include "ale3_nurbs.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 03/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Ale3Nurbs::Ale3Nurbs(int id, int owner) :
DRT::ELEMENTS::Ale3::Ale3(id,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 03/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Ale3Nurbs::Ale3Nurbs
(const DRT::ELEMENTS::NURBS::Ale3Nurbs& old) :
DRT::ELEMENTS::Ale3::Ale3(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 03/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Ale3Nurbs::~Ale3Nurbs()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 03/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::Ale3Nurbs::Print(ostream& os) const
{
  os << "Ale3Nurbs ";
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 03/09 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::NURBS::Ale3Nurbs::Shape() const
{
  switch (NumNode())
  {
  case   8: return nurbs8;
  case  27: return nurbs27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }

  return dis_none;
}

#endif // CCADISCRET

#endif //D_ALE
