/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of enrichment-wall fluid elements.
In addition to that, it contains the interface between element call
and Gauss point loop (depending on the fluid implementation)
as well as some additional service routines (for the evaluation
of errors, turbulence statistics etc.)


\level 2

*/
/*-----------------------------------------------------------*/


#include "4C_fluid_ele_xwall.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fluid_ele_nullspace.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::FluidXWallType DRT::ELEMENTS::FluidXWallType::instance_;

DRT::ELEMENTS::FluidXWallType& DRT::ELEMENTS::FluidXWallType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::FluidXWallType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::FluidXWall* object = new DRT::ELEMENTS::FluidXWall(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::FluidXWallType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDXW")
  {
    return Teuchos::rcp(new DRT::ELEMENTS::FluidXWall(id, owner));
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::FluidXWallType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::FluidXWall(id, owner));
}

void DRT::ELEMENTS::FluidXWallType::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  // this is necessary here! Otherwise it would not be consistent with the non-enriched nodes
  // since we are assuming that all elements are equal during nullspace computation
  numdf = 4;
  dimns = numdf;
  nv = numdf - 1;
  np = 1;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::FluidXWallType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

void DRT::ELEMENTS::FluidXWallType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defsxwall = definitions["FLUIDXW"];

  defsxwall["HEX8"] = INPUT::LineDefinition::Builder()
                          .AddIntVector("HEX8", 8)
                          .AddNamedInt("MAT")
                          .AddNamedString("NA")
                          .Build();
  defsxwall["TET4"] = INPUT::LineDefinition::Builder()
                          .AddIntVector("TET4", 4)
                          .AddNamedInt("MAT")
                          .AddNamedString("NA")
                          .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidXWall::FluidXWall(int id, int owner) : Fluid(id, owner) { return; }

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidXWall::FluidXWall(const DRT::ELEMENTS::FluidXWall& old) : Fluid(old) { return; }

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::FluidXWall::Clone() const
{
  DRT::ELEMENTS::FluidXWall* newelement = new DRT::ELEMENTS::FluidXWall(*this);
  return newelement;
}



/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                           |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::FluidXWall::Lines()
{
  return CORE::COMM::GetElementLines<FluidXWallBoundary, FluidXWall>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                                     |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::FluidXWall::Surfaces()
{
  return CORE::COMM::GetElementSurfaces<FluidXWallBoundary, FluidXWall>(*this);
}

/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidXWall::Print(std::ostream& os) const
{
  os << "FluidXWall ";
  Element::Print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
