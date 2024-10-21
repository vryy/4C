#include "4C_solid_poro_3D_ele_pressure_based.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"
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

namespace Discret::ELEMENTS::SolidPoroPressureBasedInternal
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
          .add_optional_named_string("TYPE");
    }
  }  // namespace
}  // namespace Discret::ELEMENTS::SolidPoroPressureBasedInternal

Discret::ELEMENTS::SolidPoroPressureBasedType
    Discret::ELEMENTS::SolidPoroPressureBasedType::instance_;

Discret::ELEMENTS::SolidPoroPressureBasedType&
Discret::ELEMENTS::SolidPoroPressureBasedType::instance()
{
  return instance_;
}

void Discret::ELEMENTS::SolidPoroPressureBasedType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsgeneral =
      definitions["SOLIDPORO_PRESSURE_BASED"];

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex8)] =
      Discret::ELEMENTS::SolidPoroPressureBasedInternal::get_default_line_definition_builder<
          Core::FE::CellType::hex8>()
          .add_optional_named_string("EAS")
          .add_optional_tag("FBAR")
          .build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex27)] =
      Discret::ELEMENTS::SolidPoroPressureBasedInternal::get_default_line_definition_builder<
          Core::FE::CellType::hex27>()
          .build();


  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::tet4)] =
      Discret::ELEMENTS::SolidPoroPressureBasedInternal::get_default_line_definition_builder<
          Core::FE::CellType::tet4>()
          .build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::tet10)] =
      Discret::ELEMENTS::SolidPoroPressureBasedInternal::get_default_line_definition_builder<
          Core::FE::CellType::tet10>()
          .build();
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidPoroPressureBasedType::create(
    const std::string eletype, const std::string elecelltype, const int id, const int owner)
{
  if (eletype == "SOLIDPORO_PRESSURE_BASED") return create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidPoroPressureBasedType::create(
    const int id, const int owner)
{
  return Teuchos::make_rcp<Discret::ELEMENTS::SolidPoroPressureBased>(id, owner);
}

Core::Communication::ParObject* Discret::ELEMENTS::SolidPoroPressureBasedType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::ELEMENTS::SolidPoroPressureBased(-1, -1);
  object->unpack(buffer);
  return object;
}

void Discret::ELEMENTS::SolidPoroPressureBasedType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  FourC::Solid::Utils::nodal_block_information_solid(dwele, numdf, dimns, nv, np);
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SolidPoroPressureBasedType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_3d_null_space(node, x0);
}

Discret::ELEMENTS::SolidPoroPressureBased::SolidPoroPressureBased(int id, int owner)
    : Core::Elements::Element(id, owner)
{
}

Core::Elements::Element* Discret::ELEMENTS::SolidPoroPressureBased::clone() const
{
  return new Discret::ELEMENTS::SolidPoroPressureBased(*this);
}

int Discret::ELEMENTS::SolidPoroPressureBased::num_line() const
{
  return Core::FE::get_number_of_element_lines(celltype_);
}

int Discret::ELEMENTS::SolidPoroPressureBased::num_surface() const
{
  return Core::FE::get_number_of_element_surfaces(celltype_);
}

int Discret::ELEMENTS::SolidPoroPressureBased::num_volume() const
{
  return Core::FE::get_number_of_element_volumes(celltype_);
}

std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::SolidPoroPressureBased::lines()
{
  return Core::Communication::get_element_lines<StructuralLine, SolidPoroPressureBased>(*this);
}

std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::SolidPoroPressureBased::surfaces()
{
  return Core::Communication::get_element_surfaces<StructuralSurface, SolidPoroPressureBased>(
      *this);
}

void Discret::ELEMENTS::SolidPoroPressureBased::set_params_interface_ptr(
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

bool Discret::ELEMENTS::SolidPoroPressureBased::read_element(const std::string& eletype,
    const std::string& elecelltype, const Core::IO::InputParameterContainer& container)
{
  // read base element
  // set cell type
  celltype_ = Core::FE::string_to_cell_type(elecelltype);

  // read number of material model
  set_material(
      0, Mat::factory(FourC::Solid::Utils::read_element::read_element_material(container)));

  // read kinematic type
  solid_ele_property_.kintype =
      FourC::Solid::Utils::read_element::read_element_kinematic_type(container);

  // check element technology
  if (FourC::Solid::Utils::read_element::read_element_technology(container) !=
      ElementTechnology::none)
    FOUR_C_THROW("SOLIDPORO elements do not support any element technology!");

  // read scalar transport implementation type
  poro_ele_property_.impltype = FourC::Solid::Utils::read_element::read_type(container);

  solid_calc_variant_ = create_solid_calculation_interface(celltype_, solid_ele_property_);
  solidporo_press_based_calc_variant_ =
      create_solid_poro_pressure_based_calculation_interface(celltype_);

  // setup solid material
  std::visit(
      [&](auto& solid) { solid->setup(struct_poro_material(), container); }, solid_calc_variant_);

  // setup poro material
  std::visit([&](auto& solidporopressurebased)
      { solidporopressurebased->poro_setup(struct_poro_material(), container); },
      solidporo_press_based_calc_variant_);

  return true;
}

Mat::So3Material& Discret::ELEMENTS::SolidPoroPressureBased::solid_poro_material(int nummat) const
{
  return *Teuchos::rcp_dynamic_cast<Mat::So3Material>(
      Core::Elements::Element::material(nummat), true);
}

void Discret::ELEMENTS::SolidPoroPressureBased::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  add_to_pack(data, unique_par_object_id());

  // add base class Element
  Core::Elements::Element::pack(data);

  add_to_pack(data, (int)celltype_);

  Discret::ELEMENTS::add_to_pack(data, solid_ele_property_);

  add_to_pack(data, poro_ele_property_.impltype);

  data.add_to_pack(material_post_setup_);

  // optional data, e.g., EAS data
  Discret::ELEMENTS::pack(solid_calc_variant_, data);
  Discret::ELEMENTS::pack(solidporo_press_based_calc_variant_, data);
}

void Discret::ELEMENTS::SolidPoroPressureBased::unpack(Core::Communication::UnpackBuffer& buffer)
{
  if (extract_int(buffer) != unique_par_object_id()) FOUR_C_THROW("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer base_buffer(basedata);
  Core::Elements::Element::unpack(base_buffer);

  celltype_ = static_cast<Core::FE::CellType>(extract_int(buffer));

  Discret::ELEMENTS::extract_from_pack(buffer, solid_ele_property_);

  poro_ele_property_.impltype = static_cast<Inpar::ScaTra::ImplType>(extract_int(buffer));

  extract_from_pack(buffer, material_post_setup_);

  // reset solid and poro interfaces
  solid_calc_variant_ = create_solid_calculation_interface(celltype_, solid_ele_property_);
  solidporo_press_based_calc_variant_ =
      create_solid_poro_pressure_based_calculation_interface(celltype_);

  Discret::ELEMENTS::unpack(solid_calc_variant_, buffer);
  Discret::ELEMENTS::unpack(solidporo_press_based_calc_variant_, buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

void Discret::ELEMENTS::SolidPoroPressureBased::vis_names(std::map<std::string, int>& names)
{
  Core::Elements::Element::vis_names(names);
  solid_poro_material().vis_names(names);
}

bool Discret::ELEMENTS::SolidPoroPressureBased::vis_data(
    const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return solid_poro_material().vis_data(name, data, id());
}

Mat::StructPoro& Discret::ELEMENTS::SolidPoroPressureBased::struct_poro_material(int nummat) const
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
Mat::FluidPoroMultiPhase& Discret::ELEMENTS::SolidPoroPressureBased::fluid_poro_material(
    int nummat) const
{
  if (this->num_material() <= 1)
  {
    FOUR_C_THROW("No second material defined for SolidPoroPressureBased element %i", id());
  }

  auto fluidmulti_mat = Teuchos::rcp_dynamic_cast<Mat::FluidPoroMultiPhase>(
      Core::Elements::Element::material(1), true);

  if (fluidmulti_mat == Teuchos::null)
    FOUR_C_THROW("cast to multiphase fluid poro material failed");
  if (fluidmulti_mat->material_type() != Core::Materials::m_fluidporo_multiphase and
      fluidmulti_mat->material_type() != Core::Materials::m_fluidporo_multiphase_reactions)
    FOUR_C_THROW("invalid fluid material for poro-multiphase-elasticity");
  if (fluidmulti_mat->num_fluid_phases() == 0)
  {
    FOUR_C_THROW(
        "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE = 0 currently not supported since this requires "
        "an adaption of the definition of the solid pressure");
  }
  return *fluidmulti_mat;
}

FOUR_C_NAMESPACE_CLOSE
