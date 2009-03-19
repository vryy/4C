/*!----------------------------------------------------------------------
\file ale2_nurbs.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#ifdef CCADISCRET

#include "ale2_nurbs.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 04/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Ale2Nurbs::Ale2Nurbs(int id, int owner) :
DRT::ELEMENTS::Ale2::Ale2(id,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Ale2Nurbs::Ale2Nurbs
(const DRT::ELEMENTS::NURBS::Ale2Nurbs& old) :
DRT::ELEMENTS::Ale2::Ale2(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 02/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Ale2Nurbs::~Ale2Nurbs()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::Ale2Nurbs::Print(ostream& os) const
{
  os << "Ale2Nurbs ";
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 02/09 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::NURBS::Ale2Nurbs::Shape() const
{
  switch (NumNode())
  {
  case  4: return nurbs4;
  case  9: return nurbs9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }

  return dis_none;
}

#endif // CCADISCRET

#endif //D_ALE
