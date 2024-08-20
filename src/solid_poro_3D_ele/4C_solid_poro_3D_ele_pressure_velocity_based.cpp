/*! \file

\brief Implementation of the solid-poro-pressure-velocity-based element

\level 1
                                         */

#include "4C_solid_poro_3D_ele_pressure_velocity_based.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_solid_3D_ele_factory.hpp"
#include "4C_solid_3D_ele_interface_serializable.hpp"
#include "4C_solid_3D_ele_utils.hpp"
#include "4C_solid_poro_3D_ele_factory.hpp"
#include "4C_solid_poro_3D_ele_utils.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS::SolidPoroPressureVelocityBasedInternal
{
  namespace
  {
    template <Core::FE::CellType celltype>
    Input::LineDefinition::Builder get_default_line_definition_builder()
    {
      return Input::LineDefinition::Builder()
          .add_int_vector(Core::FE::cell_type_to_string(celltype), Core::FE::num_nodes<celltype>)
          .add_named_int("MAT")
          .add_named_string("KINEM")
          .add_optional_named_double_vector("POROANISODIR1", 3)
          .add_optional_named_double_vector("POROANISODIR2", 3)
          .add_optional_named_double_vector("POROANISODIR3", 3);
    }
  }  // namespace
}  // namespace Discret::ELEMENTS::SolidPoroPressureVelocityBasedInternal

Discret::ELEMENTS::SolidPoroPressureVelocityBasedType
    Discret::ELEMENTS::SolidPoroPressureVelocityBasedType::instance_;

Discret::ELEMENTS::SolidPoroPressureVelocityBasedType&
Discret::ELEMENTS::SolidPoroPressureVelocityBasedType::instance()
{
  return instance_;
}

void Discret::ELEMENTS::SolidPoroPressureVelocityBasedType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsgeneral =
      definitions["SOLIDPORO_PRESSURE_VELOCITY_BASED"];

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex8)] =
      Discret::ELEMENTS::SolidPoroPressureVelocityBasedInternal::
          get_default_line_definition_builder<Core::FE::CellType::hex8>()
              .add_optional_named_double_vector("POROANISONODALCOEFFS1", 8)
              .add_optional_named_double_vector("POROANISONODALCOEFFS2", 8)
              .add_optional_named_double_vector("POROANISONODALCOEFFS3", 8)
              .build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex27)] =
      Discret::ELEMENTS::SolidPoroPressureVelocityBasedInternal::
          get_default_line_definition_builder<Core::FE::CellType::hex27>()
              .build();


  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::tet4)] =
      Discret::ELEMENTS::SolidPoroPressureVelocityBasedInternal::
          get_default_line_definition_builder<Core::FE::CellType::tet4>()
              .add_optional_named_double_vector("POROANISONODALCOEFFS1", 8)
              .add_optional_named_double_vector("POROANISONODALCOEFFS2", 8)
              .add_optional_named_double_vector("POROANISONODALCOEFFS3", 8)
              .build();


  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::tet10)] =
      Discret::ELEMENTS::SolidPoroPressureVelocityBasedInternal::
          get_default_line_definition_builder<Core::FE::CellType::tet10>()
              .build();
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidPoroPressureVelocityBasedType::create(
    const std::string eletype, const std::string elecelltype, const int id, const int owner)
{
  if (eletype == "SOLIDPORO_PRESSURE_VELOCITY_BASED") return create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidPoroPressureVelocityBasedType::create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::SolidPoroPressureVelocityBased(id, owner));
}

Core::Communication::ParObject* Discret::ELEMENTS::SolidPoroPressureVelocityBasedType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::ELEMENTS::SolidPoroPressureVelocityBased(-1, -1);
  object->unpack(buffer);
  return object;
}

void Discret::ELEMENTS::SolidPoroPressureVelocityBasedType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  FourC::Solid::UTILS::nodal_block_information_solid(dwele, numdf, dimns, nv, np);
}

Core::LinAlg::SerialDenseMatrix
Discret::ELEMENTS::SolidPoroPressureVelocityBasedType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_3d_null_space(node, x0);
}

Discret::ELEMENTS::SolidPoroPressureVelocityBased::SolidPoroPressureVelocityBased(int id, int owner)
    : Core::Elements::Element(id, owner)
{
}

Core::Elements::Element* Discret::ELEMENTS::SolidPoroPressureVelocityBased::clone() const
{
  return new Discret::ELEMENTS::SolidPoroPressureVelocityBased(*this);
}

int Discret::ELEMENTS::SolidPoroPressureVelocityBased::num_line() const
{
  return Core::FE::get_number_of_element_lines(celltype_);
}

int Discret::ELEMENTS::SolidPoroPressureVelocityBased::num_surface() const
{
  return Core::FE::get_number_of_element_surfaces(celltype_);
}

int Discret::ELEMENTS::SolidPoroPressureVelocityBased::num_volume() const
{
  return Core::FE::get_number_of_element_volumes(celltype_);
}

std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::SolidPoroPressureVelocityBased::lines()
{
  return Core::Communication::get_element_lines<StructuralLine, SolidPoroPressureVelocityBased>(
      *this);
}

std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::SolidPoroPressureVelocityBased::surfaces()
{
  return Core::Communication::get_element_surfaces<StructuralSurface,
      SolidPoroPressureVelocityBased>(*this);
}

void Discret::ELEMENTS::SolidPoroPressureVelocityBased::set_params_interface_ptr(
    const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface");
    solid_interface_ptr_ =
        Teuchos::rcp_dynamic_cast<FourC::Solid::ELEMENTS::ParamsInterface>(interface_ptr_);
  }
  else
  {
    interface_ptr_ = Teuchos::null;
    solid_interface_ptr_ = Teuchos::null;
  }
}

bool Discret::ELEMENTS::SolidPoroPressureVelocityBased::read_element(const std::string& eletype,
    const std::string& elecelltype, const Core::IO::InputParameterContainer& container)
{
  // read base element
  // set cell type
  celltype_ = Core::FE::string_to_cell_type(elecelltype);

  // set anisotropic_properties
  anisotropic_permeability_property_.directions_.resize(3);
  anisotropic_permeability_property_.nodal_coeffs_.resize(3);

  // read number of material model
  set_material(
      0, Mat::factory(FourC::Solid::UTILS::read_element::read_element_material(container)));

  // read kinematic type
  solid_ele_property_.kintype =
      FourC::Solid::UTILS::read_element::read_element_kinematic_type(container);

  // check element technology
  if (FourC::Solid::UTILS::read_element::read_element_technology(container) !=
      ElementTechnology::none)
    FOUR_C_THROW("SOLIDPORO elements do not support any element technology!");

  // read scalar transport implementation type
  poro_ele_property_.impltype = FourC::Solid::UTILS::read_element::read_type(container);

  read_anisotropic_permeability_directions_from_element_line_definition(container);
  read_anisotropic_permeability_nodal_coeffs_from_element_line_definition(container);

  solid_calc_variant_ = create_solid_calculation_interface(celltype_, solid_ele_property_);
  solidporo_press_vel_based_calc_variant_ =
      create_solid_poro_pressure_velocity_based_calculation_interface(celltype_);

  // setup solid material
  std::visit(
      [&](auto& solid) { solid->setup(struct_poro_material(), container); }, solid_calc_variant_);

  // setup poro material
  std::visit([&](auto& solidporopressurevelocitybased)
      { solidporopressurevelocitybased->poro_setup(struct_poro_material(), container); },
      solidporo_press_vel_based_calc_variant_);

  return true;
}

void Discret::ELEMENTS::SolidPoroPressureVelocityBased::
    read_anisotropic_permeability_directions_from_element_line_definition(
        const Core::IO::InputParameterContainer& container)
{
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string definition_name = "POROANISODIR" + std::to_string(dim + 1);
    anisotropic_permeability_property_.directions_[dim] =
        container.get_or<std::vector<double>>(definition_name, std::vector<double>(3, 0.0));
  }
}

void Discret::ELEMENTS::SolidPoroPressureVelocityBased::
    read_anisotropic_permeability_nodal_coeffs_from_element_line_definition(
        const Core::IO::InputParameterContainer& container)
{
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string definition_name = "POROANISONODALCOEFFS" + std::to_string(dim + 1);
    anisotropic_permeability_property_.nodal_coeffs_[dim] = container.get_or<std::vector<double>>(
        definition_name, std::vector<double>(this->num_node(), 0.0));
  }
}


Mat::So3Material& Discret::ELEMENTS::SolidPoroPressureVelocityBased::solid_poro_material(
    int nummat) const
{
  return *Teuchos::rcp_dynamic_cast<Mat::So3Material>(
      Core::Elements::Element::material(nummat), true);
}

void Discret::ELEMENTS::SolidPoroPressureVelocityBased::pack(
    Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  add_to_pack(data, unique_par_object_id());

  // add base class Element
  Core::Elements::Element::pack(data);

  add_to_pack(data, (int)celltype_);

  Discret::ELEMENTS::add_to_pack(data, solid_ele_property_);

  data.add_to_pack(material_post_setup_);

  // anisotropic_permeability_directions_
  auto size = static_cast<int>(anisotropic_permeability_property_.directions_.size());
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i)
    add_to_pack(data, anisotropic_permeability_property_.directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = static_cast<int>(anisotropic_permeability_property_.nodal_coeffs_.size());
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i)
    add_to_pack(data, anisotropic_permeability_property_.nodal_coeffs_[i]);

  // optional data, e.g., EAS data
  Discret::ELEMENTS::pack(solid_calc_variant_, data);
  Discret::ELEMENTS::pack(solidporo_press_vel_based_calc_variant_, data);
}

void Discret::ELEMENTS::SolidPoroPressureVelocityBased::unpack(
    Core::Communication::UnpackBuffer& buffer)
{
  if (extract_int(buffer) != unique_par_object_id()) FOUR_C_THROW("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer base_buffer(basedata);
  Core::Elements::Element::unpack(base_buffer);

  celltype_ = static_cast<Core::FE::CellType>(extract_int(buffer));

  Discret::ELEMENTS::extract_from_pack(buffer, solid_ele_property_);

  Core::Communication::ParObject::extract_from_pack(buffer, material_post_setup_);

  // anisotropic_permeability_directions_
  int size = 0;
  extract_from_pack(buffer, size);
  anisotropic_permeability_property_.directions_.resize(size, std::vector<double>(3, 0.0));
  for (int i = 0; i < size; ++i)
    extract_from_pack(buffer, anisotropic_permeability_property_.directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = 0;
  extract_from_pack(buffer, size);
  anisotropic_permeability_property_.nodal_coeffs_.resize(
      size, std::vector<double>(this->num_node(), 0.0));
  for (int i = 0; i < size; ++i)
    extract_from_pack(buffer, anisotropic_permeability_property_.nodal_coeffs_[i]);



  // reset solid and poro interfaces
  solid_calc_variant_ = create_solid_calculation_interface(celltype_, solid_ele_property_);
  solidporo_press_vel_based_calc_variant_ =
      create_solid_poro_pressure_velocity_based_calculation_interface(celltype_);

  Discret::ELEMENTS::unpack(solid_calc_variant_, buffer);
  Discret::ELEMENTS::unpack(solidporo_press_vel_based_calc_variant_, buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

void Discret::ELEMENTS::SolidPoroPressureVelocityBased::vis_names(std::map<std::string, int>& names)
{
  Core::Elements::Element::vis_names(names);
  solid_poro_material().vis_names(names);
}

bool Discret::ELEMENTS::SolidPoroPressureVelocityBased::vis_data(
    const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return solid_poro_material().vis_data(name, data, id());
}

Mat::StructPoro& Discret::ELEMENTS::SolidPoroPressureVelocityBased::struct_poro_material(
    int nummat) const
{
  auto porostruct_mat =
      Teuchos::rcp_dynamic_cast<Mat::StructPoro>(Core::Elements::Element::material(nummat), true);

  if (porostruct_mat == Teuchos::null) FOUR_C_THROW("cast to poro material failed");

  if (porostruct_mat->material_type() != Core::Materials::m_structporo and
      porostruct_mat->material_type() != Core::Materials::m_structpororeaction and
      porostruct_mat->material_type() != Core::Materials::m_structpororeactionECM)
    FOUR_C_THROW("invalid structure material for poroelasticity");

  return *porostruct_mat;
}


Mat::FluidPoro& Discret::ELEMENTS::SolidPoroPressureVelocityBased::fluid_poro_material(
    int nummat) const
{
  if (this->num_material() <= 1)
  {
    FOUR_C_THROW("No second material defined for SolidPoroPressureVelocityBased element %i", id());
  }

  auto fluidmulti_mat =
      Teuchos::rcp_dynamic_cast<Mat::FluidPoro>(Core::Elements::Element::material(1), true);

  if (fluidmulti_mat == Teuchos::null)
    FOUR_C_THROW("cast to multiphase fluid poro material failed");
  if (fluidmulti_mat->material_type() != Core::Materials::m_fluidporo)
    FOUR_C_THROW("invalid fluid material for poroelasticity");
  return *fluidmulti_mat;
}

FOUR_C_NAMESPACE_CLOSE
