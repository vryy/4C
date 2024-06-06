/*-----------------------------------------------------------*/
/*! \file

\brief Fluid element for poroelasticity problems


\level 2

*/
/*-----------------------------------------------------------*/


#include "4C_fluid_ele_poro.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidPoroEleType Discret::ELEMENTS::FluidPoroEleType::instance_;

Discret::ELEMENTS::FluidPoroEleType& Discret::ELEMENTS::FluidPoroEleType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::FluidPoroEleType::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::FluidPoro(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidPoroEleType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDPORO")
  {
    return Teuchos::rcp(new Discret::ELEMENTS::FluidPoro(id, owner));
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidPoroEleType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::FluidPoro(id, owner));
}

void Discret::ELEMENTS::FluidPoroEleType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_fluid;
  FluidType::setup_element_definition(definitions_fluid);

  std::map<std::string, Input::LineDefinition>& defs_fluid = definitions_fluid["FLUID"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["FLUIDPORO"];

  // 3D
  defs["HEX8"] = defs_fluid["HEX8"];
  defs["HEX20"] = defs_fluid["HEX20"];
  defs["HEX27"] = defs_fluid["HEX27"];
  defs["TET4"] = defs_fluid["TET4"];
  defs["TET10"] = defs_fluid["TET10"];
  defs["WEDGE6"] = defs_fluid["WEDGE6"];
  defs["WEDGE15"] = defs_fluid["WEDGE15"];
  defs["PYRAMID5"] = defs_fluid["PYRAMID5"];
  defs["NURBS8"] = defs_fluid["NURBS8"];
  defs["NURBS27"] = defs_fluid["NURBS27"];

  // 2D
  defs["QUAD4"] = defs_fluid["QUAD4"];
  defs["QUAD8"] = defs_fluid["QUAD8"];
  defs["QUAD9"] = defs_fluid["QUAD9"];
  defs["TRI3"] = defs_fluid["TRI3"];
  defs["TRI6"] = defs_fluid["TRI6"];
  defs["NURBS4"] = defs_fluid["NURBS4"];
  defs["NURBS9"] = defs_fluid["NURBS9"];
}

Discret::ELEMENTS::FluidPoro::FluidPoro(int id, int owner)
    : Fluid(id, owner), kintype_(Inpar::STR::KinemType::vague)
{
  anisotropic_permeability_directions_.resize(3, std::vector<double>(1, 0.0));
  anisotropic_permeability_nodal_coeffs_.resize(3, std::vector<double>(1, 0.0));
}

Discret::ELEMENTS::FluidPoro::FluidPoro(const Discret::ELEMENTS::FluidPoro& old)
    : Fluid(old),
      kintype_(old.kintype_),
      anisotropic_permeability_directions_(old.anisotropic_permeability_directions_),
      anisotropic_permeability_nodal_coeffs_(old.anisotropic_permeability_nodal_coeffs_)
{
}

Core::Elements::Element* Discret::ELEMENTS::FluidPoro::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::FluidPoro(*this);
  return newelement;
}

void Discret::ELEMENTS::FluidPoro::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // kinemtics type
  AddtoPack(data, kintype_);

  // anisotropic_permeability_directions_
  auto size = static_cast<int>(anisotropic_permeability_directions_.size());
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = static_cast<int>(anisotropic_permeability_nodal_coeffs_.size());
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, anisotropic_permeability_nodal_coeffs_[i]);

  // add base class Element
  Fluid::Pack(data);
}

void Discret::ELEMENTS::FluidPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // kintype_
  kintype_ = static_cast<Inpar::STR::KinemType>(ExtractInt(position, data));

  // anisotropic_permeability_directions_
  int size = 0;
  ExtractfromPack(position, data, size);
  anisotropic_permeability_directions_.resize(size, std::vector<double>(3, 0.0));
  for (int i = 0; i < size; ++i)
    ExtractfromPack(position, data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = 0;
  ExtractfromPack(position, data, size);
  anisotropic_permeability_nodal_coeffs_.resize(size, std::vector<double>(this->num_node(), 0.0));
  for (int i = 0; i < size; ++i)
    ExtractfromPack(position, data, anisotropic_permeability_nodal_coeffs_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  Fluid::ExtractfromPack(position, data, basedata);
  Fluid::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::FluidPoro::Lines()
{
  return Core::Communication::GetElementLines<FluidPoroBoundary, FluidPoro>(*this);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::FluidPoro::Surfaces()
{
  return Core::Communication::GetElementSurfaces<FluidPoroBoundary, FluidPoro>(*this);
}

void Discret::ELEMENTS::FluidPoro::Print(std::ostream& os) const
{
  os << "FluidPoro " << (Core::FE::CellTypeToString(distype_)).c_str();
  Element::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
