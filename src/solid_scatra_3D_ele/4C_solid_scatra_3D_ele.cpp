// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solid_scatra_3D_ele.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_solid_scatra_3D_ele_factory.hpp"
#include "4C_solid_scatra_3D_ele_lib.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <Core::FE::CellType celltype>
  Input::LineDefinition::Builder get_default_line_definition_builder()
  {
    return Input::LineDefinition::Builder()
        .add_int_vector(Core::FE::cell_type_to_string(celltype), Core::FE::num_nodes<celltype>)
        .add_named_int("MAT")
        .add_named_string("KINEM")
        .add_named_string("TYPE")
        .add_optional_named_string("PRESTRESS_TECH")
        .add_optional_named_double_vector("RAD", 3)
        .add_optional_named_double_vector("AXI", 3)
        .add_optional_named_double_vector("CIR", 3)
        .add_optional_named_double_vector("FIBER1", 3)
        .add_optional_named_double_vector("FIBER2", 3)
        .add_optional_named_double_vector("FIBER3", 3);
  }
}  // namespace

Discret::ELEMENTS::SolidScatraType Discret::ELEMENTS::SolidScatraType::instance_;

Discret::ELEMENTS::SolidScatraType& Discret::ELEMENTS::SolidScatraType::instance()
{
  return instance_;
}

void Discret::ELEMENTS::SolidScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsgeneral = definitions["SOLIDSCATRA"];

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex8)] =
      get_default_line_definition_builder<Core::FE::CellType::hex8>()
          .add_optional_named_string("TECH")
          .build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex27)] =
      get_default_line_definition_builder<Core::FE::CellType::hex27>().build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::tet4)] =
      get_default_line_definition_builder<Core::FE::CellType::tet4>().build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::tet10)] =
      get_default_line_definition_builder<Core::FE::CellType::tet10>().build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::nurbs27)] =
      get_default_line_definition_builder<Core::FE::CellType::nurbs27>().build();
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidScatraType::create(
    const std::string eletype, const std::string elecelltype, const int id, const int owner)
{
  if (eletype == "SOLIDSCATRA") return create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidScatraType::create(
    const int id, const int owner)
{
  return Teuchos::make_rcp<Discret::ELEMENTS::SolidScatra>(id, owner);
}

Core::Communication::ParObject* Discret::ELEMENTS::SolidScatraType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::ELEMENTS::SolidScatra(-1, -1);
  object->unpack(buffer);
  return object;
}

void Discret::ELEMENTS::SolidScatraType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  Solid::Utils::nodal_block_information_solid(dwele, numdf, dimns, nv, np);
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SolidScatraType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_3d_null_space(node, x0);
}

Discret::ELEMENTS::SolidScatra::SolidScatra(int id, int owner) : Core::Elements::Element(id, owner)
{
}

Core::Elements::Element* Discret::ELEMENTS::SolidScatra::clone() const
{
  return new Discret::ELEMENTS::SolidScatra(*this);
}

int Discret::ELEMENTS::SolidScatra::num_line() const
{
  return Core::FE::get_number_of_element_lines(celltype_);
}

int Discret::ELEMENTS::SolidScatra::num_surface() const
{
  return Core::FE::get_number_of_element_surfaces(celltype_);
}

int Discret::ELEMENTS::SolidScatra::num_volume() const
{
  return Core::FE::get_number_of_element_volumes(celltype_);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SolidScatra::lines()
{
  return Core::Communication::get_element_lines<StructuralLine, SolidScatra>(*this);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SolidScatra::surfaces()
{
  return Core::Communication::get_element_surfaces<StructuralSurface, SolidScatra>(*this);
}

void Discret::ELEMENTS::SolidScatra::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface");
    solid_interface_ptr_ =
        Teuchos::rcp_dynamic_cast<Solid::ELEMENTS::ParamsInterface>(interface_ptr_);
  }
  else
  {
    interface_ptr_ = Teuchos::null;
    solid_interface_ptr_ = Teuchos::null;
  }
}

bool Discret::ELEMENTS::SolidScatra::read_element(const std::string& eletype,
    const std::string& celltype, const Core::IO::InputParameterContainer& container)
{
  // read base element
  // set cell type
  celltype_ = Core::FE::string_to_cell_type(celltype);

  // read number of material model
  set_material(0, Mat::factory(Solid::Utils::read_element::read_element_material(container)));

  // read scalar transport implementation type
  properties_.impltype = read_scatra_impl_type(container);

  properties_.solid = Solid::Utils::read_element::read_solid_element_properties(container);

  solid_scatra_calc_variant_ =
      create_solid_scatra_calculation_interface(celltype_, properties_.solid);

  // setup solid material
  std::visit([&](auto& solid_scatra) { solid_scatra->setup(solid_material(), container); },
      solid_scatra_calc_variant_);

  return true;
}

void Discret::ELEMENTS::SolidScatra::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  add_to_pack(data, unique_par_object_id());

  // add base class Element
  Core::Elements::Element::pack(data);

  add_to_pack(data, (int)celltype_);
  Discret::ELEMENTS::add_to_pack(data, properties_);

  data.add_to_pack(material_post_setup_);

  // optional data, e.g., EAS data
  Discret::ELEMENTS::pack(solid_scatra_calc_variant_, data);
}

void Discret::ELEMENTS::SolidScatra::unpack(Core::Communication::UnpackBuffer& buffer)
{
  if (extract_int(buffer) != unique_par_object_id()) FOUR_C_THROW("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer base_buffer(basedata);
  Core::Elements::Element::unpack(base_buffer);

  celltype_ = static_cast<Core::FE::CellType>(extract_int(buffer));

  Discret::ELEMENTS::extract_from_pack(buffer, properties_);

  extract_from_pack(buffer, material_post_setup_);

  // reset solid and scatra interfaces
  solid_scatra_calc_variant_ =
      create_solid_scatra_calculation_interface(celltype_, properties_.solid);

  Discret::ELEMENTS::unpack(solid_scatra_calc_variant_, buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

void Discret::ELEMENTS::SolidScatra::vis_names(std::map<std::string, int>& names)
{
  Core::Elements::Element::vis_names(names);
  solid_material().vis_names(names);
}

bool Discret::ELEMENTS::SolidScatra::vis_data(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return solid_material().vis_data(name, data, id());
}

Mat::So3Material& Discret::ELEMENTS::SolidScatra::solid_material(int nummat) const
{
  return *Teuchos::rcp_dynamic_cast<Mat::So3Material>(
      Core::Elements::Element::material(nummat), true);
}

FOUR_C_NAMESPACE_CLOSE