/*! \file

\brief Implementation of the solid element

\level 1
*/


#include "4C_solid_3D_ele.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_discretization_fem_general_cell_type.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_solid_3D_ele_factory.hpp"
#include "4C_solid_3D_ele_interface_serializable.hpp"
#include "4C_solid_3D_ele_utils.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

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
        .AddOptionalNamedString("PRESTRESS_TECH")
        .AddOptionalNamedDoubleVector("RAD", 3)
        .AddOptionalNamedDoubleVector("AXI", 3)
        .AddOptionalNamedDoubleVector("CIR", 3)
        .AddOptionalNamedDoubleVector("FIBER1", 3)
        .AddOptionalNamedDoubleVector("FIBER2", 3)
        .AddOptionalNamedDoubleVector("FIBER3", 3);
  }
}  // namespace

DRT::ELEMENTS::SolidType DRT::ELEMENTS::SolidType::instance_;

DRT::ELEMENTS::SolidType& DRT::ELEMENTS::SolidType::Instance() { return instance_; }

void DRT::ELEMENTS::SolidType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defsgeneral = definitions["SOLID"];

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::hex8)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::hex8>()
          .AddOptionalNamedString("TECH")
          .Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::hex18)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::hex18>().Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::hex20)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::hex20>().Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::hex27)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::hex27>().Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::tet4)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::tet4>().Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::tet10)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::tet10>().Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::wedge6)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::wedge6>().Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::pyramid5)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::pyramid5>()
          .AddOptionalNamedString("TECH")
          .Build();



  defsgeneral["NURBS27"] = INPUT::LineDefinition::Builder()
                               .AddIntVector("NURBS27", 27)
                               .AddNamedInt("MAT")
                               .AddNamedString("KINEM")
                               .Build();
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidType::Create(
    const std::string eletype, const std::string elecelltype, const int id, const int owner)
{
  if (eletype == "SOLID") return Create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::Solid(id, owner));
}

CORE::COMM::ParObject* DRT::ELEMENTS::SolidType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Solid(-1, -1);
  object->Unpack(data);
  return object;
}

void DRT::ELEMENTS::SolidType::NodalBlockInformation(
    Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  STR::UTILS::NodalBlockInformationSolid(dwele, numdf, dimns, nv, np);
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SolidType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
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

DRT::ELEMENTS::Solid::Solid(int id, int owner) : DRT::Element(id, owner) {}


DRT::Element* DRT::ELEMENTS::Solid::Clone() const { return new Solid(*this); }

int DRT::ELEMENTS::Solid::NumLine() const { return CORE::FE::getNumberOfElementLines(celltype_); }

int DRT::ELEMENTS::Solid::NumSurface() const
{
  return CORE::FE::getNumberOfElementSurfaces(celltype_);
}

int DRT::ELEMENTS::Solid::NumVolume() const
{
  return CORE::FE::getNumberOfElementVolumes(celltype_);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Solid::Lines()
{
  return CORE::COMM::GetElementLines<StructuralLine, Solid>(*this);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Solid::Surfaces()
{
  return CORE::COMM::GetElementSurfaces<StructuralSurface, Solid>(*this);
}

void DRT::ELEMENTS::Solid::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  AddtoPack(data, UniqueParObjectId());

  // add base class Element
  DRT::Element::Pack(data);

  AddtoPack(data, (int)celltype_);

  AddToPack(data, solid_ele_property_);

  data.AddtoPack(material_post_setup_);

  DRT::ELEMENTS::Pack(solid_calc_variant_, data);
}

void DRT::ELEMENTS::Solid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  if (ExtractInt(position, data) != UniqueParObjectId()) FOUR_C_THROW("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Element::Unpack(basedata);

  celltype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));

  ExtractFromPack(position, data, solid_ele_property_);

  if (Shape() == CORE::FE::CellType::nurbs27)
  {
    SetNurbsElement() = true;
  }

  CORE::COMM::ParObject::ExtractfromPack(position, data, material_post_setup_);

  // reset solid interface
  solid_calc_variant_ = CreateSolidCalculationInterface(celltype_, solid_ele_property_);

  DRT::ELEMENTS::Unpack(solid_calc_variant_, position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

void DRT::ELEMENTS::Solid::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface"));
  }
  else
    interface_ptr_ = Teuchos::null;
}

bool DRT::ELEMENTS::Solid::ReadElement(
    const std::string& eletype, const std::string& celltype, INPUT::LineDefinition* linedef)
{
  // set cell type
  celltype_ = CORE::FE::StringToCellType(celltype);

  // read number of material model
  SetMaterial(STR::UTILS::READELEMENT::ReadElementMaterial(linedef));

  // kinematic type
  SetKinematicType(STR::UTILS::READELEMENT::ReadElementKinematicType(linedef));

  if (linedef->HaveNamed("TECH"))
  {
    solid_ele_property_.element_technology =
        STR::UTILS::READELEMENT::ReadElementTechnology(linedef);
  }

  if (linedef->HaveNamed("PRESTRESS_TECH"))
  {
    solid_ele_property_.prestress_technology =
        STR::UTILS::READELEMENT::ReadPrestressTechnology(linedef);
  }


  // kinematic type
  solid_ele_property_.kintype = STR::UTILS::READELEMENT::ReadElementKinematicType(linedef);

  if (Shape() == CORE::FE::CellType::nurbs27)
  {
    SetNurbsElement() = true;
  }

  solid_calc_variant_ = CreateSolidCalculationInterface(celltype_, solid_ele_property_);
  std::visit(
      [&](auto& interface) { interface->Setup(*SolidMaterial(), linedef); }, solid_calc_variant_);
  return true;
}

Teuchos::RCP<MAT::So3Material> DRT::ELEMENTS::Solid::SolidMaterial(int nummat) const
{
  return Teuchos::rcp_dynamic_cast<MAT::So3Material>(DRT::Element::Material(nummat), true);
}

void DRT::ELEMENTS::Solid::VisNames(std::map<std::string, int>& names)
{
  DRT::Element::VisNames(names);
  SolidMaterial()->VisNames(names);
}

bool DRT::ELEMENTS::Solid::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, Id());
}

FOUR_C_NAMESPACE_CLOSE
