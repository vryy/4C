/*! \file

\brief Implementation of the solid element

\level 1
*/

#include "baci_solid_ele.H"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_io_linedefinition.H"
#include "baci_lib_utils_factory.H"
#include "baci_mat_so3_material.H"
#include "baci_so3_line.H"
#include "baci_so3_nullspace.H"
#include "baci_so3_surface.H"
#include "baci_solid_ele_calc_interface.H"
#include "baci_solid_ele_factory.H"
#include "baci_solid_ele_interface_serializable.H"
#include "baci_solid_ele_utils.H"
#include "baci_structure_new_elements_paramsinterface.H"

#include <memory>


DRT::ELEMENTS::SolidType DRT::ELEMENTS::SolidType::instance_;

DRT::ELEMENTS::SolidType& DRT::ELEMENTS::SolidType::Instance() { return instance_; }

void DRT::ELEMENTS::SolidType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defsgeneral = definitions["SOLID"];

  defsgeneral["HEX8"] = INPUT::LineDefinition::Builder()
                            .AddIntVector("HEX8", 8)
                            .AddNamedInt("MAT")
                            .AddNamedString("KINEM")
                            .AddOptionalNamedString("EAS")
                            .AddOptionalTag("FBAR")
                            .AddOptionalTag("MULF")
                            .AddOptionalNamedDoubleVector("RAD", 3)
                            .AddOptionalNamedDoubleVector("AXI", 3)
                            .AddOptionalNamedDoubleVector("CIR", 3)
                            .AddOptionalNamedDoubleVector("FIBER1", 3)
                            .AddOptionalNamedDoubleVector("FIBER2", 3)
                            .AddOptionalNamedDoubleVector("FIBER3", 3)
                            .Build();

  defsgeneral["HEX18"] = INPUT::LineDefinition::Builder()
                             .AddIntVector("HEX18", 18)
                             .AddNamedInt("MAT")
                             .AddNamedString("KINEM")
                             .AddOptionalTag("MULF")
                             .AddOptionalNamedDoubleVector("RAD", 3)
                             .AddOptionalNamedDoubleVector("AXI", 3)
                             .AddOptionalNamedDoubleVector("CIR", 3)
                             .AddOptionalNamedDoubleVector("FIBER1", 3)
                             .AddOptionalNamedDoubleVector("FIBER2", 3)
                             .AddOptionalNamedDoubleVector("FIBER3", 3)
                             .Build();

  defsgeneral["HEX20"] = INPUT::LineDefinition::Builder()
                             .AddIntVector("HEX20", 20)
                             .AddNamedInt("MAT")
                             .AddNamedString("KINEM")
                             .AddOptionalTag("MULF")
                             .AddOptionalNamedDoubleVector("RAD", 3)
                             .AddOptionalNamedDoubleVector("AXI", 3)
                             .AddOptionalNamedDoubleVector("CIR", 3)
                             .AddOptionalNamedDoubleVector("FIBER1", 3)
                             .AddOptionalNamedDoubleVector("FIBER2", 3)
                             .AddOptionalNamedDoubleVector("FIBER3", 3)
                             .Build();

  defsgeneral["HEX27"] = INPUT::LineDefinition::Builder()
                             .AddIntVector("HEX27", 27)
                             .AddNamedInt("MAT")
                             .AddNamedString("KINEM")
                             .AddOptionalTag("MULF")
                             .AddOptionalNamedDoubleVector("RAD", 3)
                             .AddOptionalNamedDoubleVector("AXI", 3)
                             .AddOptionalNamedDoubleVector("CIR", 3)
                             .AddOptionalNamedDoubleVector("FIBER1", 3)
                             .AddOptionalNamedDoubleVector("FIBER2", 3)
                             .AddOptionalNamedDoubleVector("FIBER3", 3)
                             .Build();

  defsgeneral["NURBS27"] = INPUT::LineDefinition::Builder()
                               .AddIntVector("NURBS27", 27)
                               .AddNamedInt("MAT")
                               .AddNamedString("KINEM")
                               .Build();

  defsgeneral["TET4"] = INPUT::LineDefinition::Builder()
                            .AddIntVector("TET4", 4)
                            .AddNamedInt("MAT")
                            .AddNamedString("KINEM")
                            .AddOptionalTag("MULF")
                            .AddOptionalNamedDoubleVector("RAD", 3)
                            .AddOptionalNamedDoubleVector("AXI", 3)
                            .AddOptionalNamedDoubleVector("CIR", 3)
                            .AddOptionalNamedDoubleVector("FIBER1", 3)
                            .AddOptionalNamedDoubleVector("FIBER2", 3)
                            .AddOptionalNamedDoubleVector("FIBER3", 3)
                            .Build();

  defsgeneral["TET10"] = INPUT::LineDefinition::Builder()
                             .AddIntVector("TET10", 10)
                             .AddNamedInt("MAT")
                             .AddNamedString("KINEM")
                             .AddOptionalTag("MULF")
                             .AddOptionalNamedDoubleVector("RAD", 3)
                             .AddOptionalNamedDoubleVector("AXI", 3)
                             .AddOptionalNamedDoubleVector("CIR", 3)
                             .AddOptionalNamedDoubleVector("FIBER1", 3)
                             .AddOptionalNamedDoubleVector("FIBER2", 3)
                             .AddOptionalNamedDoubleVector("FIBER3", 3)
                             .Build();

  defsgeneral["WEDGE6"] = INPUT::LineDefinition::Builder()
                              .AddIntVector("WEDGE6", 6)
                              .AddNamedInt("MAT")
                              .AddNamedString("KINEM")
                              .AddOptionalTag("MULF")
                              .AddOptionalNamedDoubleVector("RAD", 3)
                              .AddOptionalNamedDoubleVector("AXI", 3)
                              .AddOptionalNamedDoubleVector("CIR", 3)
                              .AddOptionalNamedDoubleVector("FIBER1", 3)
                              .AddOptionalNamedDoubleVector("FIBER2", 3)
                              .AddOptionalNamedDoubleVector("FIBER3", 3)
                              .Build();

  defsgeneral["PYRAMID5"] = INPUT::LineDefinition::Builder()
                                .AddIntVector("PYRAMID5", 5)
                                .AddNamedInt("MAT")
                                .AddNamedString("KINEM")
                                .AddOptionalTag("FBAR")
                                .AddOptionalTag("MULF")
                                .AddOptionalNamedDoubleVector("RAD", 3)
                                .AddOptionalNamedDoubleVector("AXI", 3)
                                .AddOptionalNamedDoubleVector("CIR", 3)
                                .AddOptionalNamedDoubleVector("FIBER1", 3)
                                .AddOptionalNamedDoubleVector("FIBER2", 3)
                                .AddOptionalNamedDoubleVector("FIBER3", 3)
                                .Build();
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLID") return Create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::Solid(id, owner));
}

DRT::ParObject* DRT::ELEMENTS::SolidType::Create(const std::vector<char>& data)
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
      dserror(
          "The null space computation of a solid element of dimension %d is not yet implemented",
          numdof);
  }
  exit(1);
}

DRT::ELEMENTS::Solid::Solid(int id, int owner) : DRT::Element(id, owner) {}


DRT::Element* DRT::ELEMENTS::Solid::Clone() const { return new Solid(*this); }

int DRT::ELEMENTS::Solid::NumLine() const
{
  return CORE::DRT::UTILS::getNumberOfElementLines(distype_);
}

int DRT::ELEMENTS::Solid::NumSurface() const
{
  return CORE::DRT::UTILS::getNumberOfElementSurfaces(distype_);
}

int DRT::ELEMENTS::Solid::NumVolume() const
{
  return CORE::DRT::UTILS::getNumberOfElementVolumes(distype_);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Solid::Lines()
{
  return DRT::UTILS::GetElementLines<StructuralLine, Solid>(*this);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Solid::Surfaces()
{
  return DRT::UTILS::GetElementSurfaces<StructuralSurface, Solid>(*this);
}

void DRT::ELEMENTS::Solid::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  AddtoPack(data, UniqueParObjectId());

  // add base class Element
  DRT::Element::Pack(data);

  AddtoPack(data, (int)distype_);

  AddtoPack(data, (int)solid_ele_property_.kintype);

  AddtoPack(data, solid_ele_property_.eletech);

  AddtoPack(data, solid_ele_property_.eastype);

  data.AddtoPack(material_post_setup_);

  DRT::ELEMENTS::Pack(solid_calc_variant_, data);
}

void DRT::ELEMENTS::Solid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  if (ExtractInt(position, data) != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Element::Unpack(basedata);

  distype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));

  solid_ele_property_.kintype = static_cast<INPAR::STR::KinemType>(ExtractInt(position, data));

  DRT::ParObject::ExtractfromPack(position, data, solid_ele_property_.eletech);

  solid_ele_property_.eastype = static_cast<STR::ELEMENTS::EasType>(ExtractInt(position, data));

  if (Shape() == CORE::FE::CellType::nurbs27)
  {
    SetNurbsElement() = true;
  }

  DRT::ParObject::ExtractfromPack(position, data, material_post_setup_);

  // reset solid interface
  solid_calc_variant_ =
      CreateSolidCalculationInterface(*this, GetEleTech(), GetKinemType(), GetEASType());

  DRT::ELEMENTS::Unpack(solid_calc_variant_, position, data);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
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
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // set discretization type
  distype_ = CORE::FE::StringToCellType(distype);

  // read number of material model
  SetMaterial(STR::UTILS::READELEMENT::ReadElementMaterial(linedef));

  // kinematic type
  SetKinematicType(STR::UTILS::READELEMENT::ReadElementKinematicType(linedef));

  if (linedef->HaveNamed("EAS"))
  {
    if (Shape() == CORE::FE::CellType::hex8)
    {
      STR::UTILS::READELEMENT::ReadAndSetEAS(
          linedef, solid_ele_property_.eastype, solid_ele_property_.eletech);
    }
    else
      dserror("no EAS allowed for this element shape");
  }

  if (linedef->HaveNamed("FBAR"))
  {
    solid_ele_property_.eletech.insert(INPAR::STR::EleTech::fbar);
  }

  if (linedef->HaveNamed("MULF"))
  {
    solid_ele_property_.eletech.insert(INPAR::STR::EleTech::ps_mulf);
  }

  if (Shape() == CORE::FE::CellType::nurbs27)
  {
    SetNurbsElement() = true;
  }

  solid_calc_variant_ =
      CreateSolidCalculationInterface(*this, GetEleTech(), GetKinemType(), GetEASType());
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
