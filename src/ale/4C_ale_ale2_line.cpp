/*----------------------------------------------------------------------------*/
/*! \file

\brief 2D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale2.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::Ale2LineType Discret::ELEMENTS::Ale2LineType::instance_;

Discret::ELEMENTS::Ale2LineType& Discret::ELEMENTS::Ale2LineType::Instance() { return instance_; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Ale2Line::Ale2Line(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::ELEMENTS::Ale2* parent, const int lline)
    : Core::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lline);
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Ale2Line::Ale2Line(const Discret::ELEMENTS::Ale2Line& old)
    : Core::Elements::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Ale2Line::Clone() const
{
  Discret::ELEMENTS::Ale2Line* newelement = new Discret::ELEMENTS::Ale2Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Ale2Line::Shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      return Core::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale2Line::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Ale2Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale2Line::unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Ale2Line element does not support communication");
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale2Line::print(std::ostream& os) const
{
  os << "Ale2Line ";
  Element::print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
