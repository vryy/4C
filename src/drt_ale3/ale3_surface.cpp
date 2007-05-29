
#ifdef D_ALE
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "ale3.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
#include "../ale3/ale3.h"
}



DRT::Elements::Ale3Surface::Ale3Surface(int id,
                                        int owner,
                                        int nnode,
                                        const int* nodeids,
                                        DRT::Node** nodes,
                                        DRT::Elements::Ale3* parent,
                                        const int lsurface)
  : DRT::Element(id,element_ale3surface,owner),
    parent_(parent),
    lsurface_(lsurface)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
}


DRT::Elements::Ale3Surface::Ale3Surface(const DRT::Elements::Ale3Surface& old)
  : DRT::Element(old),
    parent_(old.parent_),
    lsurface_(old.lsurface_)
{
}


DRT::Element* DRT::Elements::Ale3Surface::Clone() const
{
  DRT::Elements::Ale3Surface* newelement = new DRT::Elements::Ale3Surface(*this);
  return newelement;
}


DRT::Element::DiscretizationType DRT::Elements::Ale3Surface::Shape() const
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


void DRT::Elements::Ale3Surface::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Ale3Surface element does not support communication");
}


void DRT::Elements::Ale3Surface::Unpack(const vector<char>& data)
{
  dserror("this Ale3Surface element does not support communication");
}


DRT::Elements::Ale3Surface::~Ale3Surface()
{
}


void DRT::Elements::Ale3Surface::Print(ostream& os) const
{
  os << "Ale3Surface ";
  Element::Print(os);
}


#endif
#endif
#endif
