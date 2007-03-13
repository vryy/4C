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
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3.H"
#include "linalg_utils.H"
#include "drt_utils.H"
#include "drt_discret.H"
#include "drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
}
#include "dstrc.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Surface::Fluid3Surface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::Fluid3* parent,
                              const int lsurface) :
DRT::Element(id,element_fluid3surface,owner),
parent_(parent),
lsurface_(lsurface)
{
  DSTraceHelper dst("Fluid3Surface::Fluid3Surface");
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Surface::Fluid3Surface(const DRT::Elements::Fluid3Surface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{
  DSTraceHelper dst("Fluid3Surface::Fluid3Surface");
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Fluid3Surface::Clone() const
{
  DSTraceHelper dst("Fluid3Surface::Clone");
  DRT::Elements::Fluid3Surface* newelement = new DRT::Elements::Fluid3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3Surface::Pack(vector<char>& data) const
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
void DRT::Elements::Fluid3Surface::Unpack(const vector<char>& data)
{
  DSTraceHelper dst("Fluid3Surface::Unpack");
  dserror("this Fluid3Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Surface::~Fluid3Surface()
{
  DSTraceHelper dst("Fluid3Surface::~Fluid3Surface");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3Surface::Print(ostream& os) const
{
  DSTraceHelper dst("Fluid3Surface::Print");
  os << "Fluid3Surface ";
  Element::Print(os);
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
