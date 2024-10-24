// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ale_ale3.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::Elements::Ale3SurfaceType Discret::Elements::Ale3SurfaceType::instance_;

Discret::Elements::Ale3SurfaceType& Discret::Elements::Ale3SurfaceType::instance()
{
  return instance_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::Elements::Ale3Surface::Ale3Surface(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::Elements::Ale3* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lsurface);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::Elements::Ale3Surface::Ale3Surface(const Discret::Elements::Ale3Surface& old)
    : Core::Elements::FaceElement(old)
{
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::Ale3Surface::clone() const
{
  Discret::Elements::Ale3Surface* newelement = new Discret::Elements::Ale3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::FE::CellType Discret::Elements::Ale3Surface::shape() const
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
void Discret::Elements::Ale3Surface::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Ale3Surface element does not support communication");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::Elements::Ale3Surface::unpack(Core::Communication::UnpackBuffer& buffer)
{
  FOUR_C_THROW("this Ale3Surface element does not support communication");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::Elements::Ale3Surface::print(std::ostream& os) const
{
  os << "Ale3Surface ";
  Element::print(os);
}

FOUR_C_NAMESPACE_CLOSE
