#include "4C_fluid_ele_poro.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidPoroEleType Discret::ELEMENTS::FluidPoroEleType::instance_;

Discret::ELEMENTS::FluidPoroEleType& Discret::ELEMENTS::FluidPoroEleType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::FluidPoroEleType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::ELEMENTS::FluidPoro(-1, -1);
  object->unpack(buffer);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidPoroEleType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUIDPORO")
  {
    return Teuchos::make_rcp<Discret::ELEMENTS::FluidPoro>(id, owner);
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidPoroEleType::create(
    const int id, const int owner)
{
  return Teuchos::make_rcp<Discret::ELEMENTS::FluidPoro>(id, owner);
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
    : Fluid(id, owner), kintype_(Inpar::Solid::KinemType::vague)
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

Core::Elements::Element* Discret::ELEMENTS::FluidPoro::clone() const
{
  auto* newelement = new Discret::ELEMENTS::FluidPoro(*this);
  return newelement;
}

void Discret::ELEMENTS::FluidPoro::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // kinemtics type
  add_to_pack(data, kintype_);

  // anisotropic_permeability_directions_
  auto size = static_cast<int>(anisotropic_permeability_directions_.size());
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i) add_to_pack(data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = static_cast<int>(anisotropic_permeability_nodal_coeffs_.size());
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i) add_to_pack(data, anisotropic_permeability_nodal_coeffs_[i]);

  // add base class Element
  Fluid::pack(data);
}

void Discret::ELEMENTS::FluidPoro::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // kintype_
  kintype_ = static_cast<Inpar::Solid::KinemType>(extract_int(buffer));

  // anisotropic_permeability_directions_
  int size = 0;
  extract_from_pack(buffer, size);
  anisotropic_permeability_directions_.resize(size, std::vector<double>(3, 0.0));
  for (int i = 0; i < size; ++i) extract_from_pack(buffer, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = 0;
  extract_from_pack(buffer, size);
  anisotropic_permeability_nodal_coeffs_.resize(size, std::vector<double>(this->num_node(), 0.0));
  for (int i = 0; i < size; ++i)
    extract_from_pack(buffer, anisotropic_permeability_nodal_coeffs_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  Fluid::unpack(basedata_buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::FluidPoro::lines()
{
  return Core::Communication::get_element_lines<FluidPoroBoundary, FluidPoro>(*this);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::FluidPoro::surfaces()
{
  return Core::Communication::get_element_surfaces<FluidPoroBoundary, FluidPoro>(*this);
}

void Discret::ELEMENTS::FluidPoro::print(std::ostream& os) const
{
  os << "FluidPoro " << (Core::FE::cell_type_to_string(distype_)).c_str();
  Element::print(os);
}

FOUR_C_NAMESPACE_CLOSE
