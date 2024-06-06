/*----------------------------------------------------------------------------*/
/*! \file

\brief 3D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale3.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::Ale3SurfaceType Discret::ELEMENTS::Ale3SurfaceType::instance_;

Discret::ELEMENTS::Ale3SurfaceType& Discret::ELEMENTS::Ale3SurfaceType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Ale3Surface::Ale3Surface(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::ELEMENTS::Ale3* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::ELEMENTS::Ale3Surface::Ale3Surface(const Discret::ELEMENTS::Ale3Surface& old)
    : Core::Elements::FaceElement(old)
{
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Ale3Surface::Clone() const
{
  Discret::ELEMENTS::Ale3Surface* newelement = new Discret::ELEMENTS::Ale3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Ale3Surface::Shape() const
{
  switch (num_node())
  {
    case 3:
      return Core::FE::CellType::tri3;
    case 4:
      return Core::FE::CellType::quad4;
    case 6:
      return Core::FE::CellType::tri6;
    case 8:
      return Core::FE::CellType::quad8;
    case 9:
      return Core::FE::CellType::quad9;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale3Surface::Pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Ale3Surface element does not support communication");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale3Surface::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Ale3Surface element does not support communication");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Ale3Surface::Print(std::ostream& os) const
{
  os << "Ale3Surface ";
  Element::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
