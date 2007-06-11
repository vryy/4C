/*!----------------------------------------------------------------------
\file fluid3_line.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger (Ursula)
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
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
}
#include "../drt_lib/dstrc.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Line::XFluid3Line(	int id,
                                         	int owner,
                                          int nnode,
                                          const int* nodeids,
                                         	DRT::Node** nodes,
                                         	DRT::Elements::XFluid3Surface* surfaceParent,
             										DRT::Elements::XFluid3* parent,  
                                         	const int lline) :
DRT::Element(id,element_xfluid3line,owner),
surfaceParent_(surfaceParent),
parent_(parent),
lline_(lline)
{
  DSTraceHelper dst("XFluid3Line::XFluid3Line");	
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Line::XFluid3Line(const DRT::Elements::XFluid3Line& old) :
DRT::Element(old),
surfaceParent_(old.surfaceParent_),
parent_(old.parent_),
lline_(old.lline_)
{
  DSTraceHelper dst("XFluid3Line::XFluid3Line");
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::XFluid3Line::Clone() const
{
  DSTraceHelper dst("XFluid3Line::Clone");	
  DRT::Elements::XFluid3Line* newelement = new DRT::Elements::XFluid3Line(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::XFluid3Line::Shape() const
{
  switch (NumNode())
  {
  case 2: return line2;
  case 3: return line3;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Line::Pack(vector<char>& data) const
{
  DSTraceHelper dst("XFluid3Line::Pack");
  data.resize(0);
  dserror("this Fluid3Line element does not support communication");

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Line::Unpack(const vector<char>& data)
{
  DSTraceHelper dst("XFluid3Line::Unpack");	
  dserror("this XFluid3Surface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::XFluid3Line::~XFluid3Line()
{
  DSTraceHelper dst("XFluid3Line::~XFluid3Line");	
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Line::Print(ostream& os) const
{
  DSTraceHelper dst("XFluid3Line::Print");	
  os << "XFluid3Line ";
  Element::Print(os);
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
