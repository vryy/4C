/*! \file

\brief Implementation of the solid element

\level 1
*/


#include "4C_solid_3D_ele.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_solid_3D_ele_factory.hpp"
#include "4C_solid_3D_ele_interface_serializable.hpp"
#include "4C_solid_3D_ele_properties.hpp"
#include "4C_solid_3D_ele_utils.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

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
        .add_optional_named_string("PRESTRESS_TECH")
        .add_optional_named_double_vector("RAD", 3)
        .add_optional_named_double_vector("AXI", 3)
        .add_optional_named_double_vector("CIR", 3)
        .add_optional_named_double_vector("FIBER1", 3)
        .add_optional_named_double_vector("FIBER2", 3)
        .add_optional_named_double_vector("FIBER3", 3);
  }
}  // namespace

Discret::ELEMENTS::SolidType Discret::ELEMENTS::SolidType::instance_;

Discret::ELEMENTS::SolidType& Discret::ELEMENTS::SolidType::instance() { return instance_; }

void Discret::ELEMENTS::SolidType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsgeneral = definitions["SOLID"];

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex8)] =
      get_default_line_definition_builder<Core::FE::CellType::hex8>()
          .add_optional_named_string("TECH")
          .build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex18)] =
      get_default_line_definition_builder<Core::FE::CellType::hex18>().build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex20)] =
      get_default_line_definition_builder<Core::FE::CellType::hex20>().build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::hex27)] =
      get_default_line_definition_builder<Core::FE::CellType::hex27>().build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::tet4)] =
      get_default_line_definition_builder<Core::FE::CellType::tet4>().build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::tet10)] =
      get_default_line_definition_builder<Core::FE::CellType::tet10>().build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::wedge6)] =
      get_default_line_definition_builder<Core::FE::CellType::wedge6>().build();

  defsgeneral[Core::FE::cell_type_to_string(Core::FE::CellType::pyramid5)] =
      get_default_line_definition_builder<Core::FE::CellType::pyramid5>()
          .add_optional_named_string("TECH")
          .build();



  defsgeneral["NURBS27"] = Input::LineDefinition::Builder()
                               .add_int_vector("NURBS27", 27)
                               .add_named_int("MAT")
                               .add_named_string("KINEM")
                               .build();
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidType::create(
    const std::string eletype, const std::string elecelltype, const int id, const int owner)
{
  if (eletype == "SOLID") return create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidType::create(
    const int id, const int owner)
{
  return Teuchos::make_rcp<Discret::ELEMENTS::Solid>(id, owner);
}

Core::Communication::ParObject* Discret::ELEMENTS::SolidType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::ELEMENTS::Solid(-1, -1);
  object->unpack(buffer);
  return object;
}

void Discret::ELEMENTS::SolidType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  FourC::Solid::Utils::nodal_block_information_solid(dwele, numdf, dimns, nv, np);
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SolidType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  switch (numdof)
  {
    case 3:
      return compute_solid_3d_null_space(node, x0);
    case 2:
      return compute_solid_2d_null_space(node, x0);
    default:
      FOUR_C_THROW(
          "The null space computation of a solid element of dimension %d is not yet implemented",
          numdof);
  }
  exit(1);
}

Discret::ELEMENTS::Solid::Solid(int id, int owner) : Core::Elements::Element(id, owner) {}


Core::Elements::Element* Discret::ELEMENTS::Solid::clone() const { return new Solid(*this); }

int Discret::ELEMENTS::Solid::num_line() const
{
  return Core::FE::get_number_of_element_lines(celltype_);
}

int Discret::ELEMENTS::Solid::num_surface() const
{
  return Core::FE::get_number_of_element_surfaces(celltype_);
}

int Discret::ELEMENTS::Solid::num_volume() const
{
  return Core::FE::get_number_of_element_volumes(celltype_);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Solid::lines()
{
  return Core::Communication::get_element_lines<StructuralLine, Solid>(*this);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Solid::surfaces()
{
  return Core::Communication::get_element_surfaces<StructuralSurface, Solid>(*this);
}

const Core::FE::GaussIntegration& Discret::ELEMENTS::Solid::get_gauss_rule() const
{
  return std::visit([](auto& interface) -> const Core::FE::GaussIntegration&
      { return interface->get_gauss_rule_stiffness_integration(); },
      solid_calc_variant_);
}

void Discret::ELEMENTS::Solid::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  add_to_pack(data, unique_par_object_id());

  // add base class Element
  Core::Elements::Element::pack(data);

  add_to_pack(data, (int)celltype_);

  Discret::ELEMENTS::add_to_pack(data, solid_ele_property_);

  data.add_to_pack(material_post_setup_);

  Discret::ELEMENTS::pack(solid_calc_variant_, data);
}

void Discret::ELEMENTS::Solid::unpack(Core::Communication::UnpackBuffer& buffer)
{
  if (extract_int(buffer) != unique_par_object_id()) FOUR_C_THROW("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer base_buffer(basedata);
  Core::Elements::Element::unpack(base_buffer);

  celltype_ = static_cast<Core::FE::CellType>(extract_int(buffer));

  Discret::ELEMENTS::extract_from_pack(buffer, solid_ele_property_);

  extract_from_pack(buffer, material_post_setup_);

  // reset solid interface
  solid_calc_variant_ = create_solid_calculation_interface(celltype_, solid_ele_property_);

  Discret::ELEMENTS::unpack(solid_calc_variant_, buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

void Discret::ELEMENTS::Solid::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = Teuchos::rcp_dynamic_cast<FourC::Solid::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface"));
  }
  else
    interface_ptr_ = Teuchos::null;
}

bool Discret::ELEMENTS::Solid::read_element(const std::string& eletype, const std::string& celltype,
    const Core::IO::InputParameterContainer& container)
{
  // set cell type
  celltype_ = Core::FE::string_to_cell_type(celltype);

  // read number of material model
  set_material(
      0, Mat::factory(FourC::Solid::Utils::read_element::read_element_material(container)));

  // kinematic type
  set_kinematic_type(FourC::Solid::Utils::read_element::read_element_kinematic_type(container));

  solid_ele_property_ = FourC::Solid::Utils::read_element::read_solid_element_properties(container);


  solid_calc_variant_ = create_solid_calculation_interface(celltype_, solid_ele_property_);
  std::visit([&](auto& interface) { interface->setup(*solid_material(), container); },
      solid_calc_variant_);
  return true;
}

Teuchos::RCP<Mat::So3Material> Discret::ELEMENTS::Solid::solid_material(int nummat) const
{
  return Teuchos::rcp_dynamic_cast<Mat::So3Material>(
      Core::Elements::Element::material(nummat), true);
}

void Discret::ELEMENTS::Solid::vis_names(std::map<std::string, int>& names)
{
  Core::Elements::Element::vis_names(names);
  solid_material()->vis_names(names);
}

bool Discret::ELEMENTS::Solid::vis_data(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return solid_material()->vis_data(name, data, id());
}

void Discret::ELEMENTS::Solid::for_each_gauss_point(Core::FE::Discretization& discretization,
    std::vector<int>& lm,
    const std::function<void(Mat::So3Material&, double, int)>& integrator) const
{
  std::visit(
      [&](auto& interface) {
        interface->for_each_gauss_point(*this, *solid_material(), discretization, lm, integrator);
      },
      solid_calc_variant_);
}

FOUR_C_NAMESPACE_CLOSE
