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
        .add_int_vector(Core::FE::CellTypeToString(celltype), Core::FE::num_nodes<celltype>)
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

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::hex8)] =
      get_default_line_definition_builder<Core::FE::CellType::hex8>()
          .add_optional_named_string("TECH")
          .build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::hex18)] =
      get_default_line_definition_builder<Core::FE::CellType::hex18>().build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::hex20)] =
      get_default_line_definition_builder<Core::FE::CellType::hex20>().build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::hex27)] =
      get_default_line_definition_builder<Core::FE::CellType::hex27>().build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::tet4)] =
      get_default_line_definition_builder<Core::FE::CellType::tet4>().build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::tet10)] =
      get_default_line_definition_builder<Core::FE::CellType::tet10>().build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::wedge6)] =
      get_default_line_definition_builder<Core::FE::CellType::wedge6>().build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::pyramid5)] =
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
  return Teuchos::rcp(new Discret::ELEMENTS::Solid(id, owner));
}

Core::Communication::ParObject* Discret::ELEMENTS::SolidType::create(const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::Solid(-1, -1);
  object->unpack(data);
  return object;
}

void Discret::ELEMENTS::SolidType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  FourC::Solid::UTILS::nodal_block_information_solid(dwele, numdf, dimns, nv, np);
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SolidType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  switch (numdof)
  {
    case 3:
      return ComputeSolid3DNullSpace(node, x0);
    case 2:
      return ComputeSolid2DNullSpace(node, x0);
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
  return Core::FE::getNumberOfElementLines(celltype_);
}

int Discret::ELEMENTS::Solid::num_surface() const
{
  return Core::FE::getNumberOfElementSurfaces(celltype_);
}

int Discret::ELEMENTS::Solid::num_volume() const
{
  return Core::FE::getNumberOfElementVolumes(celltype_);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Solid::lines()
{
  return Core::Communication::GetElementLines<StructuralLine, Solid>(*this);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Solid::surfaces()
{
  return Core::Communication::GetElementSurfaces<StructuralSurface, Solid>(*this);
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

void Discret::ELEMENTS::Solid::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  if (extract_int(position, data) != unique_par_object_id())
    FOUR_C_THROW("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Elements::Element::unpack(basedata);

  celltype_ = static_cast<Core::FE::CellType>(extract_int(position, data));

  Discret::ELEMENTS::ExtractFromPack(position, data, solid_ele_property_);

  if (shape() == Core::FE::CellType::nurbs27)
  {
    set_nurbs_element() = true;
  }

  Core::Communication::ParObject::extract_from_pack(position, data, material_post_setup_);

  // reset solid interface
  solid_calc_variant_ = create_solid_calculation_interface(celltype_, solid_ele_property_);

  Discret::ELEMENTS::unpack(solid_calc_variant_, position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
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
  celltype_ = Core::FE::StringToCellType(celltype);

  // read number of material model
  set_material(
      0, Mat::Factory(FourC::Solid::UTILS::read_element::read_element_material(container)));

  // kinematic type
  set_kinematic_type(FourC::Solid::UTILS::read_element::read_element_kinematic_type(container));

  solid_ele_property_ = FourC::Solid::UTILS::read_element::read_solid_element_properties(container);

  if (shape() == Core::FE::CellType::nurbs27)
  {
    set_nurbs_element() = true;
  }

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
