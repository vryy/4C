//-----------------------------------------------------------------------
/*!
\file ale3_surface.cpp

<pre>

</pre>
*/
//-----------------------------------------------------------------------
#ifdef D_ALE


#include "ale3.H"
#include "../drt_lib/drt_dserror.H"


DRT::ELEMENTS::Ale3SurfaceType DRT::ELEMENTS::Ale3SurfaceType::instance_;


DRT::ELEMENTS::Ale3Surface::Ale3Surface(int id,
                                        int owner,
                                        int nnode,
                                        const int* nodeids,
                                        DRT::Node** nodes,
                                        DRT::ELEMENTS::Ale3* parent,
                                        const int lsurface)
  : DRT::Element(id,owner),
    parent_(parent),
    lsurface_(lsurface)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
}


DRT::ELEMENTS::Ale3Surface::Ale3Surface(const DRT::ELEMENTS::Ale3Surface& old)
  : DRT::Element(old),
    parent_(old.parent_),
    lsurface_(old.lsurface_)
{
}


DRT::Element* DRT::ELEMENTS::Ale3Surface::Clone() const
{
  DRT::ELEMENTS::Ale3Surface* newelement = new DRT::ELEMENTS::Ale3Surface(*this);
  return newelement;
}


DRT::Element::DiscretizationType DRT::ELEMENTS::Ale3Surface::Shape() const
{
  switch (NumNode())
  {
  case 3: return tri3;
  case 4: return quad4;
  case 6: return tri6;
  case 8: return quad8;
  case 9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}


void DRT::ELEMENTS::Ale3Surface::Pack(DRT::PackBuffer& data) const
{
  dserror("this Ale3Surface element does not support communication");
}


void DRT::ELEMENTS::Ale3Surface::Unpack(const std::vector<char>& data)
{
  dserror("this Ale3Surface element does not support communication");
}


DRT::ELEMENTS::Ale3Surface::~Ale3Surface()
{
}


void DRT::ELEMENTS::Ale3Surface::Print(std::ostream& os) const
{
  os << "Ale3Surface ";
  Element::Print(os);
}


#endif
