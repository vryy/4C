/*----------------------------------------------------------------------------*/
/*! \file

\brief 2D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale2.hpp"
#include "4C_lib_discret.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::Ale2LineType DRT::ELEMENTS::Ale2LineType::instance_;

DRT::ELEMENTS::Ale2LineType& DRT::ELEMENTS::Ale2LineType::Instance() { return instance_; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale2Line::Ale2Line(int id, int owner, int nnode, const int* nodeids,
    CORE::Nodes::Node** nodes, DRT::ELEMENTS::Ale2* parent, const int lline)
    : CORE::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lline);
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale2Line::Ale2Line(const DRT::ELEMENTS::Ale2Line& old)
    : CORE::Elements::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::Ale2Line::Clone() const
{
  DRT::ELEMENTS::Ale2Line* newelement = new DRT::ELEMENTS::Ale2Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Ale2Line::Shape() const
{
  switch (num_node())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2Line::Pack(CORE::COMM::PackBuffer& data) const
{
  FOUR_C_THROW("this Ale2Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2Line::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Ale2Line element does not support communication");
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2Line::Print(std::ostream& os) const
{
  os << "Ale2Line ";
  Element::Print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
