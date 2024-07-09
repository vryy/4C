/*! \file

\brief Implementation of the solid-poro element

\level 1
*/

#include "4C_solid_poro_3D_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
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

namespace
{
  template <Core::FE::CellType celltype>
  Input::LineDefinition::Builder GetDefaultLineDefinitionBuilder()
  {
    return Input::LineDefinition::Builder()
        .add_int_vector(Core::FE::CellTypeToString(celltype), Core::FE::num_nodes<celltype>)
        .add_named_int("MAT")
        .add_named_string("KINEM")
        .add_optional_named_double_vector("RAD", 3)
        .add_optional_named_double_vector("AXI", 3)
        .add_optional_named_double_vector("CIR", 3)
        .add_optional_named_double_vector("FIBER1", 3)
        .add_optional_named_double_vector("FIBER2", 3)
        .add_optional_named_double_vector("FIBER3", 3)
        .add_optional_named_string("TYPE")
        .add_optional_named_string("POROTYPE");
  }
}  // namespace

Discret::ELEMENTS::SolidPoroType Discret::ELEMENTS::SolidPoroType::instance_;

Discret::ELEMENTS::SolidPoroType& Discret::ELEMENTS::SolidPoroType::instance() { return instance_; }

void Discret::ELEMENTS::SolidPoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsgeneral = definitions["SOLIDPORO"];

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::hex8)] =
      GetDefaultLineDefinitionBuilder<Core::FE::CellType::hex8>()
          .add_optional_named_string("EAS")
          .add_optional_tag("FBAR")
          .build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::hex27)] =
      GetDefaultLineDefinitionBuilder<Core::FE::CellType::hex27>().build();


  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::tet4)] =
      GetDefaultLineDefinitionBuilder<Core::FE::CellType::tet4>().build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::tet10)] =
      GetDefaultLineDefinitionBuilder<Core::FE::CellType::tet10>().build();
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidPoroType::create(
    const std::string eletype, const std::string elecelltype, const int id, const int owner)
{
  if (eletype == "SOLIDPORO") return create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidPoroType::create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::SolidPoro(id, owner));
}

Core::Communication::ParObject* Discret::ELEMENTS::SolidPoroType::create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::SolidPoro(-1, -1);
  object->unpack(data);
  return object;
}

void Discret::ELEMENTS::SolidPoroType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  FourC::Solid::UTILS::nodal_block_information_solid(dwele, numdf, dimns, nv, np);
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SolidPoroType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

Discret::ELEMENTS::SolidPoro::SolidPoro(int id, int owner) : Core::Elements::Element(id, owner) {}

Core::Elements::Element* Discret::ELEMENTS::SolidPoro::clone() const
{
  return new Discret::ELEMENTS::SolidPoro(*this);
}

int Discret::ELEMENTS::SolidPoro::num_line() const
{
  return Core::FE::getNumberOfElementLines(celltype_);
}

int Discret::ELEMENTS::SolidPoro::num_surface() const
{
  return Core::FE::getNumberOfElementSurfaces(celltype_);
}

int Discret::ELEMENTS::SolidPoro::num_volume() const
{
  return Core::FE::getNumberOfElementVolumes(celltype_);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SolidPoro::lines()
{
  return Core::Communication::GetElementLines<StructuralLine, SolidPoro>(*this);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SolidPoro::surfaces()
{
  return Core::Communication::GetElementSurfaces<StructuralSurface, SolidPoro>(*this);
}

void Discret::ELEMENTS::SolidPoro::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = Teuchos::rcp_dynamic_cast<FourC::Solid::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface"));
  }
  else
    interface_ptr_ = Teuchos::null;
}

bool Discret::ELEMENTS::SolidPoro::read_element(
    const std::string& eletype, const std::string& elecelltype, Input::LineDefinition* linedef)
{
  // read base element
  // set cell type
  celltype_ = Core::FE::StringToCellType(elecelltype);

  // read number of material model
  set_material(0, Mat::Factory(FourC::Solid::UTILS::ReadElement::read_element_material(linedef)));

  // kinematic type
  solid_ele_property_.kintype =
      FourC::Solid::UTILS::ReadElement::read_element_kinematic_type(linedef);

  // check element technology
  if (linedef->has_named("TECH"))
  {
    if (FourC::Solid::UTILS::ReadElement::read_element_technology(linedef) !=
        ElementTechnology::none)
      FOUR_C_THROW("SOLIDPORO elements do not support any element technology!");
  }

  // read scalar transport implementation type
  if (linedef->has_named("POROTYPE"))
  {
    poro_ele_property_.porotype = FourC::Solid::UTILS::ReadElement::ReadPoroType(linedef);
  }
  else
  {
    poro_ele_property_.porotype = Inpar::Poro::PoroType::undefined;
  }

  // read scalar transport implementation type
  if (linedef->has_named("TYPE"))
  {
    poro_ele_property_.impltype = FourC::Solid::UTILS::ReadElement::read_type(linedef);
  }
  else
  {
    poro_ele_property_.impltype = Inpar::ScaTra::impltype_undefined;
  }

  solid_calc_variant_ = create_solid_calculation_interface(celltype_, solid_ele_property_);
  solidporo_calc_variant_ = create_solid_poro_calculation_interface(*this, get_ele_poro_type());

  // setup solid material
  std::visit(
      [&](auto& solid) { solid->setup(struct_poro_material(), linedef); }, solid_calc_variant_);

  // setup poro material
  std::visit([&](auto& solidporo) { solidporo->poro_setup(struct_poro_material(), linedef); },
      solidporo_calc_variant_);

  return true;
}

Mat::So3Material& Discret::ELEMENTS::SolidPoro::solid_poro_material(int nummat) const
{
  return *Teuchos::rcp_dynamic_cast<Mat::So3Material>(
      Core::Elements::Element::material(nummat), true);
}

void Discret::ELEMENTS::SolidPoro::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  add_to_pack(data, unique_par_object_id());

  // add base class Element
  Core::Elements::Element::pack(data);

  add_to_pack(data, (int)celltype_);

  Discret::ELEMENTS::add_to_pack(data, solid_ele_property_);

  add_to_pack(data, poro_ele_property_.porotype);

  add_to_pack(data, poro_ele_property_.impltype);

  data.add_to_pack(material_post_setup_);

  // optional data, e.g., EAS data
  Discret::ELEMENTS::pack(solid_calc_variant_, data);
  Discret::ELEMENTS::pack(solidporo_calc_variant_, data);
}

void Discret::ELEMENTS::SolidPoro::unpack(const std::vector<char>& data)
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

  poro_ele_property_.porotype = static_cast<Inpar::Poro::PoroType>(extract_int(position, data));

  poro_ele_property_.impltype = static_cast<Inpar::ScaTra::ImplType>(extract_int(position, data));

  Core::Communication::ParObject::extract_from_pack(position, data, material_post_setup_);

  // reset solid and poro interfaces
  solid_calc_variant_ = create_solid_calculation_interface(celltype_, solid_ele_property_);
  solidporo_calc_variant_ = create_solid_poro_calculation_interface(*this, get_ele_poro_type());

  Discret::ELEMENTS::unpack(solid_calc_variant_, position, data);
  Discret::ELEMENTS::unpack(solidporo_calc_variant_, position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

void Discret::ELEMENTS::SolidPoro::vis_names(std::map<std::string, int>& names)
{
  Core::Elements::Element::vis_names(names);
  solid_poro_material().vis_names(names);
}

bool Discret::ELEMENTS::SolidPoro::vis_data(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return solid_poro_material().vis_data(name, data, id());
}

Mat::StructPoro& Discret::ELEMENTS::SolidPoro::struct_poro_material(int nummat) const
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
Mat::FluidPoroMultiPhase& Discret::ELEMENTS::SolidPoro::fluid_poro_multi_material(int nummat) const
{
  if (this->num_material() <= 1)
  {
    FOUR_C_THROW("No second material defined for SolidPoro element %i", id());
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
