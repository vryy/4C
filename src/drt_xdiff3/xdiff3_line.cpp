/*!----------------------------------------------------------------------
\file xfluid3_line.cpp
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

#include "xdiff3.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XDiff3Line::XDiff3Line(
    const int id,
    const int owner,
    const int nnode,
    const int* nodeids,
    DRT::Node** nodes,
    DRT::Element* parent,
    const int lline) :
DRT::Element(id,element_xdiff3line,owner),
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
DRT::ELEMENTS::XDiff3Line::XDiff3Line(const DRT::ELEMENTS::XDiff3Line& old) :
DRT::Element(old),
parent_(old.parent_),
lline_(old.lline_)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::XDiff3Line::Clone() const
{
  DRT::ELEMENTS::XDiff3Line* newelement = new DRT::ELEMENTS::XDiff3Line(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::XDiff3Line::Shape() const
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
void DRT::ELEMENTS::XDiff3Line::Pack(std::vector<char>& data) const
{
  data.resize(0);
  dserror("this XDiff3Line element does not support communication");

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XDiff3Line::Unpack(const std::vector<char>& data)
{
  dserror("this XDiff3Line element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XDiff3Line::~XDiff3Line()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::XDiff3Line::Print(ostream& os) const
{
  os << "XDiff3Line ";
  Element::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
