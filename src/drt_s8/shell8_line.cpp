/*----------------------------------------------------------------------------*/
/*! \file
\brief shell8

\level 1

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/
#include "shell8.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"


DRT::ELEMENTS::Shell8LineType DRT::ELEMENTS::Shell8LineType::instance_;

DRT::ELEMENTS::Shell8LineType& DRT::ELEMENTS::Shell8LineType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8Line::Shell8Line(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Shell8* parent, const int lline)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lline);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8Line::Shell8Line(const DRT::ELEMENTS::Shell8Line& old) : DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Shell8Line::Clone() const
{
  DRT::ELEMENTS::Shell8Line* newelement = new DRT::ELEMENTS::Shell8Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Shell8Line::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return line2;
    case 3:
      return line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8Line::Pack(DRT::PackBuffer& data) const
{
  dserror("this Shell8Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8Line::Unpack(const std::vector<char>& data)
{
  dserror("this line element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8Line::~Shell8Line() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8Line::Print(std::ostream& os) const
{
  os << "Shell8Line ";
  Element::Print(os);
  return;
}
