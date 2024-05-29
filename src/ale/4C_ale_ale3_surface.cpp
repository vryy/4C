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


DRT::ELEMENTS::Ale3SurfaceType DRT::ELEMENTS::Ale3SurfaceType::instance_;

DRT::ELEMENTS::Ale3SurfaceType& DRT::ELEMENTS::Ale3SurfaceType::Instance() { return instance_; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale3Surface::Ale3Surface(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Ale3* parent, const int lsurface)
    : CORE::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
DRT::ELEMENTS::Ale3Surface::Ale3Surface(const DRT::ELEMENTS::Ale3Surface& old)
    : CORE::Elements::FaceElement(old)
{
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::Ale3Surface::Clone() const
{
  DRT::ELEMENTS::Ale3Surface* newelement = new DRT::ELEMENTS::Ale3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Ale3Surface::Shape() const
{
  switch (num_node())
  {
    case 3:
      return CORE::FE::CellType::tri3;
    case 4:
      return CORE::FE::CellType::quad4;
    case 6:
      return CORE::FE::CellType::tri6;
    case 8:
      return CORE::FE::CellType::quad8;
    case 9:
      return CORE::FE::CellType::quad9;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale3Surface::Pack(CORE::COMM::PackBuffer& data) const
{
  FOUR_C_THROW("this Ale3Surface element does not support communication");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale3Surface::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Ale3Surface element does not support communication");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale3Surface::Print(std::ostream& os) const
{
  os << "Ale3Surface ";
  Element::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
