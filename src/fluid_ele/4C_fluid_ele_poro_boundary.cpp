// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_ele_poro.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::FluidPoroBoundaryType Discret::Elements::FluidPoroBoundaryType::instance_;

Discret::Elements::FluidPoroBoundaryType& Discret::Elements::FluidPoroBoundaryType::instance()
{
  return instance_;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::FluidPoroBoundaryType::create(
    const int id, const int owner)
{
  return Teuchos::null;
}

Discret::Elements::FluidPoroBoundary::FluidPoroBoundary(int id, int owner, int nnode,
    const int* nodeids, Core::Nodes::Node** nodes, Discret::Elements::Fluid* parent,
    const int lsurface)
    : FluidBoundary(id, owner, nnode, nodeids, nodes, parent, lsurface)
{
}

Discret::Elements::FluidPoroBoundary::FluidPoroBoundary(
    const Discret::Elements::FluidPoroBoundary& old)
    : FluidBoundary(old)
{
}

Core::Elements::Element* Discret::Elements::FluidPoroBoundary::clone() const
{
  auto* newelement = new Discret::Elements::FluidPoroBoundary(*this);
  return newelement;
}

void Discret::Elements::FluidPoroBoundary::print(std::ostream& os) const
{
  os << "FluidPoroBoundary ";
  Element::print(os);
}

FOUR_C_NAMESPACE_CLOSE
