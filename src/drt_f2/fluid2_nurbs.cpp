/*!----------------------------------------------------------------------
\file fluid2_nurbs.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "fluid2_nurbs.H"

DRT::ELEMENTS::NURBS::Fluid2NurbsType DRT::ELEMENTS::NURBS::Fluid2NurbsType::instance_;


DRT::ParObject* DRT::ELEMENTS::NURBS::Fluid2NurbsType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::NURBS::Fluid2Nurbs* object = new DRT::ELEMENTS::NURBS::Fluid2Nurbs(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::Fluid2NurbsType::Create( const string eletype,
                                                                          const string eledistype,
                                                                          const int id,
                                                                          const int owner )
{
  if ( eletype=="FLUID2" )
  {
    if ( eledistype=="NURBS4" || eledistype=="NURBS9" )
    {
      return rcp(new DRT::ELEMENTS::NURBS::Fluid2Nurbs(id,owner));
    }
  }
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Fluid2Nurbs::Fluid2Nurbs(int id, int owner) :
DRT::ELEMENTS::Fluid2::Fluid2(id,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Fluid2Nurbs::Fluid2Nurbs
(const DRT::ELEMENTS::NURBS::Fluid2Nurbs& old) :
DRT::ELEMENTS::Fluid2::Fluid2(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::Fluid2Nurbs::~Fluid2Nurbs()
{
  return;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          gammi 05/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::NURBS::Fluid2Nurbs::Shape() const
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

#endif

#endif //D_FLUID2
