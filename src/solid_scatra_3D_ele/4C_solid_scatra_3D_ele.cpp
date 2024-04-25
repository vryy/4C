/*! \file

\brief Implementation of the solid-scatra element

\level 1
*/

#include "4C_solid_scatra_3D_ele.hpp"

#include "4C_discretization_fem_general_cell_type.hpp"
#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_solid_scatra_3D_ele_lib.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <CORE::FE::CellType celltype>
  INPUT::LineDefinition::Builder GetDefaultLineDefinitionBuilder()
  {
    return INPUT::LineDefinition::Builder()
        .AddIntVector(CORE::FE::CellTypeToString(celltype), CORE::FE::num_nodes<celltype>)
        .AddNamedInt("MAT")
        .AddNamedString("KINEM")
        .AddNamedString("TYPE")
        .AddOptionalNamedDoubleVector("RAD", 3)
        .AddOptionalNamedDoubleVector("AXI", 3)
        .AddOptionalNamedDoubleVector("CIR", 3)
        .AddOptionalNamedDoubleVector("FIBER1", 3)
        .AddOptionalNamedDoubleVector("FIBER2", 3)
        .AddOptionalNamedDoubleVector("FIBER3", 3);
  }
}  // namespace

DRT::ELEMENTS::SolidScatraType DRT::ELEMENTS::SolidScatraType::instance_;

DRT::ELEMENTS::SolidScatraType& DRT::ELEMENTS::SolidScatraType::Instance() { return instance_; }

void DRT::ELEMENTS::SolidScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defsgeneral = definitions["SOLIDSCATRA"];

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::hex8)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::hex8>().Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::hex27)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::hex27>().Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::tet4)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::tet4>().Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::tet10)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::tet10>().Build();
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidScatraType::Create(
    const std::string eletype, const std::string elecelltype, const int id, const int owner)
{
  if (eletype == "SOLIDSCATRA") return Create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidScatraType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::SolidScatra(id, owner));
}

CORE::COMM::ParObject* DRT::ELEMENTS::SolidScatraType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SolidScatra(-1, -1);
  object->Unpack(data);
  return object;
}

void DRT::ELEMENTS::SolidScatraType::NodalBlockInformation(
    Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  STR::UTILS::NodalBlockInformationSolid(dwele, numdf, dimns, nv, np);
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SolidScatraType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

DRT::ELEMENTS::SolidScatra::SolidScatra(int id, int owner) : DRT::Element(id, owner) {}

DRT::Element* DRT::ELEMENTS::SolidScatra::Clone() const
{
  return new DRT::ELEMENTS::SolidScatra(*this);
}

int DRT::ELEMENTS::SolidScatra::NumLine() const
{
  return CORE::FE::getNumberOfElementLines(celltype_);
}

int DRT::ELEMENTS::SolidScatra::NumSurface() const
{
  return CORE::FE::getNumberOfElementSurfaces(celltype_);
}

int DRT::ELEMENTS::SolidScatra::NumVolume() const
{
  return CORE::FE::getNumberOfElementVolumes(celltype_);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::SolidScatra::Lines()
{
  return CORE::COMM::GetElementLines<StructuralLine, SolidScatra>(*this);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::SolidScatra::Surfaces()
{
  return CORE::COMM::GetElementSurfaces<StructuralSurface, SolidScatra>(*this);
}

void DRT::ELEMENTS::SolidScatra::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface"));
  }
  else
    interface_ptr_ = Teuchos::null;
}

bool DRT::ELEMENTS::SolidScatra::ReadElement(
    const std::string& eletype, const std::string& celltype, INPUT::LineDefinition* linedef)
{
  // read base element
  // set cell type
  celltype_ = CORE::FE::StringToCellType(celltype);

  // read number of material model
  SetMaterial(STR::UTILS::READELEMENT::ReadElementMaterial(linedef));

  // read scalar transport implementation type
  properties_.impltype = ReadScatraImplType(*linedef);

  solid_scatra_calc_variant_ = CreateSolidScatraCalculationInterface(celltype_);

  // setup solid material
  std::visit([&](auto& solid_scatra) { solid_scatra->Setup(SolidMaterial(), linedef); },
      solid_scatra_calc_variant_);

  return true;
}

void DRT::ELEMENTS::SolidScatra::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  AddtoPack(data, UniqueParObjectId());

  // add base class Element
  DRT::Element::Pack(data);

  AddtoPack(data, (int)celltype_);

  AddtoPack(data, properties_.impltype);

  data.AddtoPack(material_post_setup_);

  // optional data, e.g., EAS data
  DRT::ELEMENTS::Pack(solid_scatra_calc_variant_, data);
}

void DRT::ELEMENTS::SolidScatra::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  if (ExtractInt(position, data) != UniqueParObjectId()) FOUR_C_THROW("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Element::Unpack(basedata);

  celltype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));

  properties_.impltype = static_cast<INPAR::SCATRA::ImplType>(ExtractInt(position, data));

  CORE::COMM::ParObject::ExtractfromPack(position, data, material_post_setup_);

  // reset solid and scatra interfaces
  solid_scatra_calc_variant_ = CreateSolidScatraCalculationInterface(celltype_);

  DRT::ELEMENTS::Unpack(solid_scatra_calc_variant_, position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

void DRT::ELEMENTS::SolidScatra::VisNames(std::map<std::string, int>& names)
{
  DRT::Element::VisNames(names);
  SolidMaterial().VisNames(names);
}

bool DRT::ELEMENTS::SolidScatra::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return SolidMaterial().VisData(name, data, Id());
}

MAT::So3Material& DRT::ELEMENTS::SolidScatra::SolidMaterial(int nummat) const
{
  return *Teuchos::rcp_dynamic_cast<MAT::So3Material>(DRT::Element::Material(nummat), true);
}

FOUR_C_NAMESPACE_CLOSE