/*!----------------------------------------------------------------------
\file combust3_line.cpp
\brief

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#include "combust3.H"


DRT::ELEMENTS::Combust3LineType DRT::ELEMENTS::Combust3LineType::instance_;

DRT::ELEMENTS::Combust3LineType& DRT::ELEMENTS::Combust3LineType::Instance()
{
  return instance_;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3Line::Combust3Line(
    const int id,
    const int owner,
    const int nnode,
    const int* nodeids,
    DRT::Node** nodes,
    DRT::Element* parent,
    const int lline) :
DRT::FaceElement(id,owner)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent,lline);
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3Line::Combust3Line(const DRT::ELEMENTS::Combust3Line& old) :
DRT::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Combust3Line::Clone() const
{
  DRT::ELEMENTS::Combust3Line* newelement = new DRT::ELEMENTS::Combust3Line(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Combust3Line::Shape() const
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
void DRT::ELEMENTS::Combust3Line::Pack(std::vector<char>& data) const
{
  dserror("this Combust3Line element does not support communication");

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3Line::Unpack(const std::vector<char>& data)
{
  dserror("this Combust3Line element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3Line::~Combust3Line()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Combust3Line::Print(std::ostream& os) const
{
  os << "Combust3Line ";
  Element::Print(os);
  return;
}

