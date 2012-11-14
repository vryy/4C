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

#include "ale2_nurbs.H"

DRT::ELEMENTS::NURBS::Ale2_NurbsType DRT::ELEMENTS::NURBS::Ale2_NurbsType::instance_;


DRT::ParObject* DRT::ELEMENTS::NURBS::Ale2_NurbsType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::NURBS::Ale2Nurbs* object = new DRT::ELEMENTS::NURBS::Ale2Nurbs(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::Ale2_NurbsType::Create( const string eletype,
                                                                         const string eledistype,
                                                                         const int id,
                                                                         const int owner )
{
  if ( eletype=="ALE2" )
  {
    if(eledistype=="NURBS4" || eledistype=="NURBS9")
    {
      return Teuchos::rcp(new DRT::ELEMENTS::NURBS::Ale2Nurbs(id,owner));
    }
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::Ale2_NurbsType::Create( const int id, const int owner )
{
  return Teuchos::rcp(new DRT::ELEMENTS::NURBS::Ale2Nurbs(id,owner));
}


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


#endif //D_ALE
