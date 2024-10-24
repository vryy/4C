// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_shell7p_line.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::Shell7pLineType Discret::Elements::Shell7pLineType::instance_;

Discret::Elements::Shell7pLineType& Discret::Elements::Shell7pLineType::instance()
{
  return instance_;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::Shell7pLineType::create(
    const int id, const int owner)
{
  return Teuchos::null;
}

Discret::Elements::Shell7pLine::Shell7pLine(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Core::Elements::Element* parent, const int lline)
    : Core::Elements::FaceElement(id, owner)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lline);
  // type of gaussian integration
  switch (shape())
  {
    case Core::FE::CellType::line2:
      gaussrule_ = Core::FE::GaussRule1D::line_2point;
      break;
    case Core::FE::CellType::line3:
      gaussrule_ = Core::FE::GaussRule1D::line_3point;
      break;
    default:
      FOUR_C_THROW("shape type unknown!\n");
  }
}

Discret::Elements::Shell7pLine::Shell7pLine(const Discret::Elements::Shell7pLine& old)
    : Core::Elements::FaceElement(old), gaussrule_(old.gaussrule_)
{
}

Core::Elements::Element* Discret::Elements::Shell7pLine::clone() const
{
  auto* newelement = new Discret::Elements::Shell7pLine(*this);
  return newelement;
}

Core::FE::CellType Discret::Elements::Shell7pLine::shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      return Core::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}

void Discret::Elements::Shell7pLine::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Shell7line element does not support communication");
}


void Discret::Elements::Shell7pLine::unpack(Core::Communication::UnpackBuffer& buffer)
{
  FOUR_C_THROW("this Shell line element does not support communication");
}


void Discret::Elements::Shell7pLine::print(std::ostream& os) const
{
  os << "Shell7pLine ";
  Element::print(os);
}

FOUR_C_NAMESPACE_CLOSE
