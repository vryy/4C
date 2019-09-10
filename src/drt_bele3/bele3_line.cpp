/*----------------------------------------------------------------------*/
/*! \file

\brief dummy 3D boundary element without any physics

\maintainer Amadeus Gebauer

\level 2
*/
/*----------------------------------------------------------------------*/

#include "bele3.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"


DRT::ELEMENTS::Bele3LineType DRT::ELEMENTS::Bele3LineType::instance_;

DRT::ELEMENTS::Bele3LineType& DRT::ELEMENTS::Bele3LineType::Instance() { return instance_; }


/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3Line::Bele3Line(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Bele3* parent, const int lline)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lline);
  SetNumDofPerNode(parent->NumDofPerNode(*nodes[0]));
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3Line::Bele3Line(const DRT::ELEMENTS::Bele3Line& old)
    : DRT::FaceElement(old), numdofpernode_(old.numdofpernode_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Bele3Line::Clone() const
{
  DRT::ELEMENTS::Bele3Line* newelement = new DRT::ELEMENTS::Bele3Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Bele3Line::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return line2;
    case 3:
      return line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
      break;
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3Line::Pack(DRT::PackBuffer& data) const
{
  dserror("this Bele3Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3Line::Unpack(const std::vector<char>& data)
{
  dserror("this Bele3Line element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3Line::~Bele3Line() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3Line::Print(std::ostream& os) const
{
  os << "Bele3_" << numdofpernode_ << "Line ";
  Element::Print(os);
  return;
}
