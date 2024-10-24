// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ale_ale2.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::Elements::Ale2LineType Discret::Elements::Ale2LineType::instance_;

Discret::Elements::Ale2LineType& Discret::Elements::Ale2LineType::instance() { return instance_; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::Elements::Ale2Line::Ale2Line(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::Elements::Ale2* parent, const int lline)
    : Core::Elements::FaceElement(id, owner)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lline);
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Discret::Elements::Ale2Line::Ale2Line(const Discret::Elements::Ale2Line& old)
    : Core::Elements::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::Ale2Line::clone() const
{
  Discret::Elements::Ale2Line* newelement = new Discret::Elements::Ale2Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Core::FE::CellType Discret::Elements::Ale2Line::shape() const
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
void Discret::Elements::Ale2Line::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Ale2Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::Elements::Ale2Line::unpack(Core::Communication::UnpackBuffer& buffer)
{
  FOUR_C_THROW("this Ale2Line element does not support communication");
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Discret::Elements::Ale2Line::print(std::ostream& os) const
{
  os << "Ale2Line ";
  Element::print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
