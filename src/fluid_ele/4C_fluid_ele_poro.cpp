// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_ele_poro.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::FluidPoroEleType Discret::Elements::FluidPoroEleType::instance_;

Discret::Elements::FluidPoroEleType& Discret::Elements::FluidPoroEleType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::FluidPoroEleType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::FluidPoro(-1, -1);
  object->unpack(buffer);
  return object;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::FluidPoroEleType::create(
    const std::string& eletype, Core::FE::CellType celltype, const int id, const int owner)
{
  if (eletype == "FLUIDPORO")
  {
    return std::make_shared<Discret::Elements::FluidPoro>(id, owner);
  }
  return nullptr;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::FluidPoroEleType::create(
    const int id, const int owner)
{
  return std::make_shared<Discret::Elements::FluidPoro>(id, owner);
}

void Discret::Elements::FluidPoroEleType::setup_element_definition(
    std::map<std::string, std::map<Core::FE::CellType, Core::IO::InputSpec>>& definitions)
{
  std::map<std::string, std::map<Core::FE::CellType, Core::IO::InputSpec>> definitions_fluid;
  FluidType::setup_element_definition(definitions_fluid);

  auto& defs_fluid = definitions_fluid["FLUID"];

  auto& defs = definitions["FLUIDPORO"];

  using namespace Core::IO::InputSpecBuilders;

  // 3D
  defs[Core::FE::CellType::hex8] = defs_fluid[Core::FE::CellType::hex8];
  defs[Core::FE::CellType::hex20] = defs_fluid[Core::FE::CellType::hex20];
  defs[Core::FE::CellType::hex27] = defs_fluid[Core::FE::CellType::hex27];
  defs[Core::FE::CellType::tet4] = defs_fluid[Core::FE::CellType::tet4];
  defs[Core::FE::CellType::tet10] = defs_fluid[Core::FE::CellType::tet10];
  defs[Core::FE::CellType::wedge6] = defs_fluid[Core::FE::CellType::wedge6];
  defs[Core::FE::CellType::wedge15] = defs_fluid[Core::FE::CellType::wedge15];
  defs[Core::FE::CellType::pyramid5] = defs_fluid[Core::FE::CellType::pyramid5];
  defs[Core::FE::CellType::nurbs8] = defs_fluid[Core::FE::CellType::nurbs8];
  defs[Core::FE::CellType::nurbs27] = defs_fluid[Core::FE::CellType::nurbs27];

  // 2D
  defs[Core::FE::CellType::quad4] = defs_fluid[Core::FE::CellType::quad4];
  defs[Core::FE::CellType::quad8] = defs_fluid[Core::FE::CellType::quad8];
  defs[Core::FE::CellType::quad9] = defs_fluid[Core::FE::CellType::quad9];
  defs[Core::FE::CellType::tri3] = defs_fluid[Core::FE::CellType::tri3];
  defs[Core::FE::CellType::tri6] = defs_fluid[Core::FE::CellType::tri6];
  defs[Core::FE::CellType::nurbs4] = defs_fluid[Core::FE::CellType::nurbs4];
  defs[Core::FE::CellType::nurbs9] = defs_fluid[Core::FE::CellType::nurbs9];
}

Discret::Elements::FluidPoro::FluidPoro(int id, int owner)
    : Fluid(id, owner), kintype_(Inpar::Solid::KinemType::vague)
{
  anisotropic_permeability_directions_.resize(3, std::vector<double>(1, 0.0));
  anisotropic_permeability_nodal_coeffs_.resize(3, std::vector<double>(1, 0.0));
}

Discret::Elements::FluidPoro::FluidPoro(const Discret::Elements::FluidPoro& old)
    : Fluid(old),
      kintype_(old.kintype_),
      anisotropic_permeability_directions_(old.anisotropic_permeability_directions_),
      anisotropic_permeability_nodal_coeffs_(old.anisotropic_permeability_nodal_coeffs_)
{
}

Core::Elements::Element* Discret::Elements::FluidPoro::clone() const
{
  auto* newelement = new Discret::Elements::FluidPoro(*this);
  return newelement;
}

void Discret::Elements::FluidPoro::pack(Core::Communication::PackBuffer& data) const
{
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

void Discret::Elements::FluidPoro::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // kintype_
  extract_from_pack(buffer, kintype_);

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
  Fluid::unpack(buffer);
}

std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::FluidPoro::lines()
{
  return Core::Communication::get_element_lines<FluidPoroBoundary, FluidPoro>(*this);
}

std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::FluidPoro::surfaces()
{
  return Core::Communication::get_element_surfaces<FluidPoroBoundary, FluidPoro>(*this);
}

void Discret::Elements::FluidPoro::print(std::ostream& os) const
{
  os << "FluidPoro " << (Core::FE::cell_type_to_string(distype_)).c_str();
  Element::print(os);
}

FOUR_C_NAMESPACE_CLOSE
