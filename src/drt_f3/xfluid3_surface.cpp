/*!----------------------------------------------------------------------
\file xfluid3_surface.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "xfluid3.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3Surface::XFluid3Surface(
    int id,
    int owner,
    const int nnode,
    const int* nodeids,
    DRT::Node** nodes,
    DRT::ELEMENTS::XFluid3* parent,
    const int lsurface) :
  DRT::Element(id, element_xfluid3surface, owner), parent_(parent),
      lsurface_(lsurface)
{
  lines_.resize(0);
  lineptrs_.resize(0);
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3Surface::XFluid3Surface(const DRT::ELEMENTS::XFluid3Surface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::XFluid3Surface::Clone() const
{
  DRT::ELEMENTS::XFluid3Surface* newelement = new DRT::ELEMENTS::XFluid3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::XFluid3Surface::Shape() const
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

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3Surface::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this XFluid3Surface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3Surface::Unpack(const vector<char>& data)
{
  dserror("this XFluid3Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XFluid3Surface::~XFluid3Surface()
{
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::XFluid3Surface::Lines()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)lines_.size()    == NumLine() &&
      (int)lineptrs_.size() == NumLine() &&
      dynamic_cast<DRT::ELEMENTS::XFluid3Line*>(lineptrs_[0]) )
    return (DRT::Element**)(&(lineptrs_[0]));
  
  // so we have to allocate new line elements
  DRT::UTILS::ElementBoundaryFactory<XFluid3Line,XFluid3Surface>(DRT::UTILS::buildLines,lines_,lineptrs_,this);

  return (DRT::Element**)(&(lineptrs_[0]));

}
  

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XFluid3Surface::Print(ostream& os) const
{
  os << "XFluid3Surface ";
  Element::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
