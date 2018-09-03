/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_poro_boundary.cpp

 \brief  the poro boundary element

\level 2

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
 *----------------------------------------------------------------------*/

#include "fluid_ele_poro.H"

DRT::ELEMENTS::FluidPoroBoundaryType DRT::ELEMENTS::FluidPoroBoundaryType::instance_;

DRT::ELEMENTS::FluidPoroBoundaryType& DRT::ELEMENTS::FluidPoroBoundaryType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::FluidPoroBoundaryType::Create(
    const int id, const int owner)
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidPoroBoundary::FluidPoroBoundary(int id, int owner, int nnode,
    const int* nodeids, DRT::Node** nodes, DRT::ELEMENTS::Fluid* parent, const int lsurface)
    : FluidBoundary(id, owner, nnode, nodeids, nodes, parent, lsurface)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidPoroBoundary::FluidPoroBoundary(const DRT::ELEMENTS::FluidPoroBoundary& old)
    : FluidBoundary(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::FluidPoroBoundary::Clone() const
{
  DRT::ELEMENTS::FluidPoroBoundary* newelement = new DRT::ELEMENTS::FluidPoroBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoroBoundary::Print(std::ostream& os) const
{
  os << "FluidPoroBoundary ";
  Element::Print(os);
  return;
}
