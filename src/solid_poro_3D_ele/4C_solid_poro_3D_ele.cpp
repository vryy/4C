/*! \file

\brief Implementation of the solid-poro element

\level 1
*/

#include "4C_solid_poro_3D_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
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
        .AddIntVector(Core::FE::CellTypeToString(celltype), Core::FE::num_nodes<celltype>)
        .AddNamedInt("MAT")
        .AddNamedString("KINEM")
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

Discret::ELEMENTS::SolidPoroType& Discret::ELEMENTS::SolidPoroType::Instance() { return instance_; }

void Discret::ELEMENTS::SolidPoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsgeneral = definitions["SOLIDPORO"];

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::hex8)] =
      GetDefaultLineDefinitionBuilder<Core::FE::CellType::hex8>()
          .add_optional_named_string("EAS")
          .AddOptionalTag("FBAR")
          .Build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::hex27)] =
      GetDefaultLineDefinitionBuilder<Core::FE::CellType::hex27>().Build();


  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::tet4)] =
      GetDefaultLineDefinitionBuilder<Core::FE::CellType::tet4>().Build();

  defsgeneral[Core::FE::CellTypeToString(Core::FE::CellType::tet10)] =
      GetDefaultLineDefinitionBuilder<Core::FE::CellType::tet10>().Build();
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidPoroType::Create(
    const std::string eletype, const std::string elecelltype, const int id, const int owner)
{
  if (eletype == "SOLIDPORO") return Create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SolidPoroType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::SolidPoro(id, owner));
}

Core::Communication::ParObject* Discret::ELEMENTS::SolidPoroType::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::SolidPoro(-1, -1);
  object->Unpack(data);
  return object;
}

void Discret::ELEMENTS::SolidPoroType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  STR::UTILS::NodalBlockInformationSolid(dwele, numdf, dimns, nv, np);
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SolidPoroType::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

Discret::ELEMENTS::SolidPoro::SolidPoro(int id, int owner) : Core::Elements::Element(id, owner) {}

Core::Elements::Element* Discret::ELEMENTS::SolidPoro::Clone() const
{
  return new Discret::ELEMENTS::SolidPoro(*this);
}

int Discret::ELEMENTS::SolidPoro::NumLine() const
{
  return Core::FE::getNumberOfElementLines(celltype_);
}

int Discret::ELEMENTS::SolidPoro::NumSurface() const
{
  return Core::FE::getNumberOfElementSurfaces(celltype_);
}

int Discret::ELEMENTS::SolidPoro::NumVolume() const
{
  return Core::FE::getNumberOfElementVolumes(celltype_);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SolidPoro::Lines()
{
  return Core::Communication::GetElementLines<StructuralLine, SolidPoro>(*this);
}

std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SolidPoro::Surfaces()
{
  return Core::Communication::GetElementSurfaces<StructuralSurface, SolidPoro>(*this);
}

void Discret::ELEMENTS::SolidPoro::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface"));
  }
  else
    interface_ptr_ = Teuchos::null;
}

bool Discret::ELEMENTS::SolidPoro::ReadElement(
    const std::string& eletype, const std::string& elecelltype, Input::LineDefinition* linedef)
{
  // read base element
  // set cell type
  celltype_ = Core::FE::StringToCellType(elecelltype);

  // read number of material model
  SetMaterial(0, Mat::Factory(STR::UTILS::ReadElement::ReadElementMaterial(linedef)));

  // kinematic type
  solid_ele_property_.kintype = STR::UTILS::ReadElement::ReadElementKinematicType(linedef);

  // check element technology
  if (linedef->HaveNamed("TECH"))
  {
    if (STR::UTILS::ReadElement::ReadElementTechnology(linedef) != ElementTechnology::none)
      FOUR_C_THROW("SOLIDPORO elements do not support any element technology!");
  }

  // read scalar transport implementation type
  if (linedef->HaveNamed("POROTYPE"))
  {
    poro_ele_property_.porotype = STR::UTILS::ReadElement::ReadPoroType(linedef);
  }
  else
  {
    poro_ele_property_.porotype = Inpar::Poro::PoroType::undefined;
  }

  // read scalar transport implementation type
  if (linedef->HaveNamed("TYPE"))
  {
    poro_ele_property_.impltype = STR::UTILS::ReadElement::ReadType(linedef);
  }
  else
  {
    poro_ele_property_.impltype = Inpar::ScaTra::impltype_undefined;
  }

  solid_calc_variant_ = CreateSolidCalculationInterface(celltype_, solid_ele_property_);
  solidporo_calc_variant_ = CreateSolidPoroCalculationInterface(*this, GetElePoroType());

  // setup solid material
  std::visit(
      [&](auto& solid) { solid->Setup(StructPoroMaterial(), linedef); }, solid_calc_variant_);

  // setup poro material
  std::visit([&](auto& solidporo) { solidporo->PoroSetup(StructPoroMaterial(), linedef); },
      solidporo_calc_variant_);

  return true;
}

Mat::So3Material& Discret::ELEMENTS::SolidPoro::SolidPoroMaterial(int nummat) const
{
  return *Teuchos::rcp_dynamic_cast<Mat::So3Material>(
      Core::Elements::Element::Material(nummat), true);
}

void Discret::ELEMENTS::SolidPoro::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  add_to_pack(data, UniqueParObjectId());

  // add base class Element
  Core::Elements::Element::Pack(data);

  add_to_pack(data, (int)celltype_);

  AddToPack(data, solid_ele_property_);

  add_to_pack(data, poro_ele_property_.porotype);

  add_to_pack(data, poro_ele_property_.impltype);

  data.add_to_pack(material_post_setup_);

  // optional data, e.g., EAS data
  Discret::ELEMENTS::Pack(solid_calc_variant_, data);
  Discret::ELEMENTS::Pack(solidporo_calc_variant_, data);
}

void Discret::ELEMENTS::SolidPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  if (ExtractInt(position, data) != UniqueParObjectId()) FOUR_C_THROW("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Elements::Element::Unpack(basedata);

  celltype_ = static_cast<Core::FE::CellType>(ExtractInt(position, data));

  ExtractFromPack(position, data, solid_ele_property_);

  poro_ele_property_.porotype = static_cast<Inpar::Poro::PoroType>(ExtractInt(position, data));

  poro_ele_property_.impltype = static_cast<Inpar::ScaTra::ImplType>(ExtractInt(position, data));

  Core::Communication::ParObject::extract_from_pack(position, data, material_post_setup_);

  // reset solid and poro interfaces
  solid_calc_variant_ = CreateSolidCalculationInterface(celltype_, solid_ele_property_);
  solidporo_calc_variant_ = CreateSolidPoroCalculationInterface(*this, GetElePoroType());

  Discret::ELEMENTS::Unpack(solid_calc_variant_, position, data);
  Discret::ELEMENTS::Unpack(solidporo_calc_variant_, position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

void Discret::ELEMENTS::SolidPoro::VisNames(std::map<std::string, int>& names)
{
  Core::Elements::Element::VisNames(names);
  SolidPoroMaterial().VisNames(names);
}

bool Discret::ELEMENTS::SolidPoro::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::VisData(name, data)) return true;

  return SolidPoroMaterial().VisData(name, data, Id());
}

Mat::StructPoro& Discret::ELEMENTS::SolidPoro::StructPoroMaterial(int nummat) const
{
  auto porostruct_mat =
      Teuchos::rcp_dynamic_cast<Mat::StructPoro>(Core::Elements::Element::Material(nummat), true);

  if (porostruct_mat == Teuchos::null) FOUR_C_THROW("cast to poro material failed");

  if (porostruct_mat->MaterialType() != Core::Materials::m_structporo and
      porostruct_mat->MaterialType() != Core::Materials::m_structpororeaction and
      porostruct_mat->MaterialType() != Core::Materials::m_structpororeactionECM)
    FOUR_C_THROW("invalid structure material for poroelasticity");

  return *porostruct_mat;
}
Mat::FluidPoroMultiPhase& Discret::ELEMENTS::SolidPoro::fluid_poro_multi_material(int nummat) const
{
  if (this->NumMaterial() <= 1)
  {
    FOUR_C_THROW("No second material defined for SolidPoro element %i", Id());
  }

  auto fluidmulti_mat = Teuchos::rcp_dynamic_cast<Mat::FluidPoroMultiPhase>(
      Core::Elements::Element::Material(1), true);

  if (fluidmulti_mat == Teuchos::null)
    FOUR_C_THROW("cast to multiphase fluid poro material failed");
  if (fluidmulti_mat->MaterialType() != Core::Materials::m_fluidporo_multiphase and
      fluidmulti_mat->MaterialType() != Core::Materials::m_fluidporo_multiphase_reactions)
    FOUR_C_THROW("invalid fluid material for poro-multiphase-elasticity");
  if (fluidmulti_mat->NumFluidPhases() == 0)
  {
    FOUR_C_THROW(
        "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE = 0 currently not supported since this requires "
        "an adaption of the definition of the solid pressure");
  }
  return *fluidmulti_mat;
}

FOUR_C_NAMESPACE_CLOSE
