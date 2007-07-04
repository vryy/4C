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
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Line::Fluid3Line(  int id,
                                        int owner,
                                        int nnode,
                                        const int* nodeids,
                                        DRT::Node** nodes,
                                        DRT::Elements::Fluid3Surface* surfaceParent,
                                        DRT::Elements::Fluid3* parent,  
                                        const int lline) :
DRT::Element(id,element_fluid3line,owner),
surfaceParent_(surfaceParent),
parent_(parent),
lline_(lline)
{
    SetNodeIds(nnode,nodeids);
    BuildNodalPointers(nodes);
    return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Line::Fluid3Line(const DRT::Elements::Fluid3Line& old) :
DRT::Element(old),
surfaceParent_(old.surfaceParent_),
parent_(old.parent_),
lline_(old.lline_)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Fluid3Line::Clone() const
{
  DRT::Elements::Fluid3Line* newelement = new DRT::Elements::Fluid3Line(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Fluid3Line::Shape() const
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
void DRT::Elements::Fluid3Line::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Fluid3Line element does not support communication");

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3Line::Unpack(const vector<char>& data)
{
  dserror("this Fluid3Line element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Fluid3Line::~Fluid3Line()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Fluid3Line::Print(ostream& os) const
{
  os << "Fluid3Line ";
  Element::Print(os);
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
