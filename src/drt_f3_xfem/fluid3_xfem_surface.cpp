/*!----------------------------------------------------------------------
\file fluid3_surface.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3_xfem.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
}
#include "../drt_lib/dstrc.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Surface::XFluid3Surface(int id,
                                              int owner,
                                              int nnode,
                                              const int* nodeids,
                                              DRT::Node** nodes,
                                              DRT::Elements::XFluid3* parent,
                                              const int lsurface) :
DRT::Element(id,element_xfluid3surface,owner),
parent_(parent),
lsurface_(lsurface)
{
  DSTraceHelper dst("XFluid3Surface::XFluid3Surface");
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Surface::XFluid3Surface(const DRT::Elements::XFluid3Surface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{
  DSTraceHelper dst("XFluid3Surface::XFluid3Surface");
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::XFluid3Surface::Clone() const
{
  DSTraceHelper dst("XFluid3Surface::Clone");
  DRT::Elements::XFluid3Surface* newelement = new DRT::Elements::XFluid3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::XFluid3Surface::Shape() const
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
void DRT::Elements::XFluid3Surface::Pack(vector<char>& data) const
{
  DSTraceHelper dst("Fluid3Surface::Pack");
  data.resize(0);
  dserror("this Fluid3Surface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Surface::Unpack(const vector<char>& data)
{
  DSTraceHelper dst("XFluid3Surface::Unpack");
  dserror("this XFluid3Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Surface::~XFluid3Surface()
{
  DSTraceHelper dst("XFluid3Surface::~XFluid3Surface");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Surface::Print(ostream& os) const
{
  DSTraceHelper dst("XFluid3Surface::Print");
  os << "XFluid3Surface ";
  Element::Print(os);
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
