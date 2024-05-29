/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of a boundary element for fluid poro problems


\level 2

*/
/*-----------------------------------------------------------*/


#include "4C_fluid_ele_poro.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::FluidPoroBoundaryType DRT::ELEMENTS::FluidPoroBoundaryType::instance_;

DRT::ELEMENTS::FluidPoroBoundaryType& DRT::ELEMENTS::FluidPoroBoundaryType::Instance()
{
  return instance_;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::FluidPoroBoundaryType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}

DRT::ELEMENTS::FluidPoroBoundary::FluidPoroBoundary(int id, int owner, int nnode,
    const int* nodeids, CORE::Nodes::Node** nodes, DRT::ELEMENTS::Fluid* parent, const int lsurface)
    : FluidBoundary(id, owner, nnode, nodeids, nodes, parent, lsurface)
{
}

DRT::ELEMENTS::FluidPoroBoundary::FluidPoroBoundary(const DRT::ELEMENTS::FluidPoroBoundary& old)
    : FluidBoundary(old)
{
}

CORE::Elements::Element* DRT::ELEMENTS::FluidPoroBoundary::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::FluidPoroBoundary(*this);
  return newelement;
}

void DRT::ELEMENTS::FluidPoroBoundary::Print(std::ostream& os) const
{
  os << "FluidPoroBoundary ";
  Element::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
