/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "vele3.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"



DRT::ELEMENTS::Vele3LineType DRT::ELEMENTS::Vele3LineType::instance_;

DRT::ELEMENTS::Vele3LineType& DRT::ELEMENTS::Vele3LineType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Line::Vele3Line(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::Element* parent, const int lline)
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
DRT::ELEMENTS::Vele3Line::Vele3Line(const DRT::ELEMENTS::Vele3Line& old) : DRT::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Vele3Line::Clone() const
{
  DRT::ELEMENTS::Vele3Line* newelement = new DRT::ELEMENTS::Vele3Line(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Vele3Line::Shape() const
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
void DRT::ELEMENTS::Vele3Line::Pack(DRT::PackBuffer& data) const
{
  dserror("this Vele3Line element does not support communication");

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Line::Unpack(const std::vector<char>& data)
{
  dserror("this Vele3Line element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Line::~Vele3Line() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Line::Print(std::ostream& os) const
{
  os << "Vele3Line ";
  Element::Print(os);
  return;
}
