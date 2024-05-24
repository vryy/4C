/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_bele_vele3.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN



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
  set_parent_master_element(parent, lline);
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
CORE::FE::CellType DRT::ELEMENTS::Vele3Line::Shape() const
{
  switch (num_node())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Line::Pack(CORE::COMM::PackBuffer& data) const
{
  FOUR_C_THROW("this Vele3Line element does not support communication");

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Line::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Vele3Line element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Line::Print(std::ostream& os) const
{
  os << "Vele3Line ";
  Element::Print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
