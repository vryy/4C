/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of a boundary element for fluid poro problems


\level 2

*/
/*-----------------------------------------------------------*/


#include "baci_fluid_ele_poro.H"

BACI_NAMESPACE_OPEN

DRT::ELEMENTS::FluidPoroBoundaryType DRT::ELEMENTS::FluidPoroBoundaryType::instance_;

DRT::ELEMENTS::FluidPoroBoundaryType& DRT::ELEMENTS::FluidPoroBoundaryType::Instance()
{
  return instance_;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidPoroBoundaryType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}

DRT::ELEMENTS::FluidPoroBoundary::FluidPoroBoundary(int id, int owner, int nnode,
    const int* nodeids, DRT::Node** nodes, DRT::ELEMENTS::Fluid* parent, const int lsurface)
    : FluidBoundary(id, owner, nnode, nodeids, nodes, parent, lsurface)
{
}

DRT::ELEMENTS::FluidPoroBoundary::FluidPoroBoundary(const DRT::ELEMENTS::FluidPoroBoundary& old)
    : FluidBoundary(old)
{
}

DRT::Element* DRT::ELEMENTS::FluidPoroBoundary::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::FluidPoroBoundary(*this);
  return newelement;
}

void DRT::ELEMENTS::FluidPoroBoundary::Print(std::ostream& os) const
{
  os << "FluidPoroBoundary ";
  Element::Print(os);
}

BACI_NAMESPACE_CLOSE
