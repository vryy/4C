/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of a boundary element for fluid poro problems


\level 2

*/
/*-----------------------------------------------------------*/


#include "4C_fluid_ele_poro.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidPoroBoundaryType Discret::ELEMENTS::FluidPoroBoundaryType::instance_;

Discret::ELEMENTS::FluidPoroBoundaryType& Discret::ELEMENTS::FluidPoroBoundaryType::Instance()
{
  return instance_;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidPoroBoundaryType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}

Discret::ELEMENTS::FluidPoroBoundary::FluidPoroBoundary(int id, int owner, int nnode,
    const int* nodeids, Core::Nodes::Node** nodes, Discret::ELEMENTS::Fluid* parent,
    const int lsurface)
    : FluidBoundary(id, owner, nnode, nodeids, nodes, parent, lsurface)
{
}

Discret::ELEMENTS::FluidPoroBoundary::FluidPoroBoundary(
    const Discret::ELEMENTS::FluidPoroBoundary& old)
    : FluidBoundary(old)
{
}

Core::Elements::Element* Discret::ELEMENTS::FluidPoroBoundary::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::FluidPoroBoundary(*this);
  return newelement;
}

void Discret::ELEMENTS::FluidPoroBoundary::Print(std::ostream& os) const
{
  os << "FluidPoroBoundary ";
  Element::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
