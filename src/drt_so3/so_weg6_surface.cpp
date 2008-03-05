/*!----------------------------------------------------------------------
\file so_weg6_surface.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_weg6.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sow6Surface::Sow6Surface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::So_weg6* parent,
                              const int lsurface) :
DRT::Element(id,element_sow6surface,owner),
parent_(parent),
lsurface_(lsurface)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sow6Surface::Sow6Surface(const DRT::ELEMENTS::Sow6Surface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Sow6Surface::Clone() const
{
  DRT::ELEMENTS::Sow6Surface* newelement = new DRT::ELEMENTS::Sow6Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Sow6Surface::Shape() const
{
  const int numnode = this->NumNode();
  switch (numnode){
  case 3: return tri3;
  break;
  case 4: return quad4;
  break;
  default:  dserror("Can not assign surface shape of wedge element");
  return quad4;
  }
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sow6Surface::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Sow6Surface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sow6Surface::Unpack(const vector<char>& data)
{
  dserror("this Sow6Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sow6Surface::~Sow6Surface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sow6Surface::Print(ostream& os) const
{
  os << "Sow6Surface ";
  Element::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOLID3
