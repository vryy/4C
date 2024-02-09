/*! \file

\brief Implementation of the solid-poro element

\level 1
*/

#include "baci_solid_poro_ele.hpp"

#include "baci_comm_utils_factory.hpp"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_mat_fluidporo_multiphase.hpp"
#include "baci_mat_structporo.hpp"
#include "baci_so3_line.hpp"
#include "baci_so3_nullspace.hpp"
#include "baci_so3_surface.hpp"
#include "baci_solid_ele_factory.hpp"
#include "baci_solid_ele_interface_serializable.hpp"
#include "baci_solid_ele_utils.hpp"
#include "baci_solid_poro_ele_factory.hpp"
#include "baci_solid_poro_ele_utils.hpp"

#include <memory>

BACI_NAMESPACE_OPEN

namespace
{
  template <CORE::FE::CellType celltype>
  INPUT::LineDefinition::Builder GetDefaultLineDefinitionBuilder()
  {
    return INPUT::LineDefinition::Builder()
        .AddIntVector(CORE::FE::CellTypeToString(celltype), CORE::FE::num_nodes<celltype>)
        .AddNamedInt("MAT")
        .AddNamedString("KINEM")
        .AddOptionalNamedDoubleVector("RAD", 3)
        .AddOptionalNamedDoubleVector("AXI", 3)
        .AddOptionalNamedDoubleVector("CIR", 3)
        .AddOptionalNamedDoubleVector("FIBER1", 3)
        .AddOptionalNamedDoubleVector("FIBER2", 3)
        .AddOptionalNamedDoubleVector("FIBER3", 3)
        .AddOptionalNamedString("TYPE")
        .AddOptionalNamedString("POROTYPE")
        .AddOptionalNamedDoubleVector("POROANISODIR1", 3)
        .AddOptionalNamedDoubleVector("POROANISODIR2", 3)
        .AddOptionalNamedDoubleVector("POROANISODIR3", 3);
  }
}  // namespace

DRT::ELEMENTS::SolidPoroType DRT::ELEMENTS::SolidPoroType::instance_;

DRT::ELEMENTS::SolidPoroType& DRT::ELEMENTS::SolidPoroType::Instance() { return instance_; }

void DRT::ELEMENTS::SolidPoroType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defsgeneral = definitions["SOLIDPORO"];

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::hex8)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::hex8>()
          .AddOptionalNamedString("EAS")
          .AddOptionalTag("FBAR")
          .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS1", 8)
          .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS2", 8)
          .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS3", 8)
          .Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::hex27)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::hex27>().Build();


  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::tet4)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::tet4>()
          .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS1", 8)
          .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS2", 8)
          .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS3", 8)
          .Build();

  defsgeneral[CORE::FE::CellTypeToString(CORE::FE::CellType::tet10)] =
      GetDefaultLineDefinitionBuilder<CORE::FE::CellType::tet10>().Build();
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidPoroType::Create(
    const std::string eletype, const std::string elecelltype, const int id, const int owner)
{
  if (eletype == "SOLIDPORO") return Create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidPoroType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::SolidPoro(id, owner));
}

CORE::COMM::ParObject* DRT::ELEMENTS::SolidPoroType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SolidPoro(-1, -1);
  object->Unpack(data);
  return object;
}

void DRT::ELEMENTS::SolidPoroType::NodalBlockInformation(
    Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  STR::UTILS::NodalBlockInformationSolid(dwele, numdf, dimns, nv, np);
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SolidPoroType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

DRT::ELEMENTS::SolidPoro::SolidPoro(int id, int owner) : DRT::Element(id, owner) {}

DRT::Element* DRT::ELEMENTS::SolidPoro::Clone() const
{
  return new DRT::ELEMENTS::SolidPoro(*this);
}

int DRT::ELEMENTS::SolidPoro::NumLine() const
{
  return CORE::FE::getNumberOfElementLines(celltype_);
}

int DRT::ELEMENTS::SolidPoro::NumSurface() const
{
  return CORE::FE::getNumberOfElementSurfaces(celltype_);
}

int DRT::ELEMENTS::SolidPoro::NumVolume() const
{
  return CORE::FE::getNumberOfElementVolumes(celltype_);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::SolidPoro::Lines()
{
  return CORE::COMM::GetElementLines<StructuralLine, SolidPoro>(*this);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::SolidPoro::Surfaces()
{
  return CORE::COMM::GetElementSurfaces<StructuralSurface, SolidPoro>(*this);
}

void DRT::ELEMENTS::SolidPoro::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface"));
  }
  else
    interface_ptr_ = Teuchos::null;
}

bool DRT::ELEMENTS::SolidPoro::ReadElement(
    const std::string& eletype, const std::string& elecelltype, INPUT::LineDefinition* linedef)
{
  // read base element
  // set cell type
  celltype_ = CORE::FE::StringToCellType(elecelltype);
  anisotropic_permeability_directions_.resize(3, std::vector<double>(3, 0.0));
  anisotropic_permeability_nodal_coeffs_.resize(3, std::vector<double>(this->NumNode(), 0.0));

  // read number of material model
  SetMaterial(STR::UTILS::READELEMENT::ReadElementMaterial(linedef));

  // kinematic type
  solid_ele_property_.kintype = STR::UTILS::READELEMENT::ReadElementKinematicType(linedef);


  if (linedef->HaveNamed("TECH"))
  {
    solid_ele_property_.element_technology =
        STR::UTILS::READELEMENT::ReadElementTechnology(linedef);
  }

  // read scalar transport implementation type
  if (linedef->HaveNamed("POROTYPE"))
  {
    poro_ele_property_.porotype = STR::UTILS::READELEMENT::ReadPoroType(linedef);
  }
  else
  {
    poro_ele_property_.porotype = INPAR::PORO::PoroType::undefined;
  }

  // read scalar transport implementation type
  if (linedef->HaveNamed("TYPE"))
  {
    poro_ele_property_.impltype = STR::UTILS::READELEMENT::ReadType(linedef);
  }
  else
  {
    poro_ele_property_.impltype = INPAR::SCATRA::impltype_undefined;
  }

  ReadAnisotropicPermeabilityDirectionsFromElementLineDefinition(linedef);
  ReadAnisotropicPermeabilityNodalCoeffsFromElementLineDefinition(linedef);

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

void DRT::ELEMENTS::SolidPoro::ReadAnisotropicPermeabilityDirectionsFromElementLineDefinition(
    INPUT::LineDefinition* linedef)
{
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string definition_name = "POROANISODIR" + std::to_string(dim + 1);
    if (linedef->HaveNamed(definition_name))
      linedef->ExtractDoubleVector(definition_name, anisotropic_permeability_directions_[dim]);
  }
}

void DRT::ELEMENTS::SolidPoro::ReadAnisotropicPermeabilityNodalCoeffsFromElementLineDefinition(
    INPUT::LineDefinition* linedef)
{
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string definition_name = "POROANISONODALCOEFFS" + std::to_string(dim + 1);
    if (linedef->HaveNamed(definition_name))
      linedef->ExtractDoubleVector(definition_name, anisotropic_permeability_nodal_coeffs_[dim]);
  }
}

MAT::So3Material& DRT::ELEMENTS::SolidPoro::SolidPoroMaterial(int nummat) const
{
  return *Teuchos::rcp_dynamic_cast<MAT::So3Material>(DRT::Element::Material(nummat), true);
}

void DRT::ELEMENTS::SolidPoro::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  AddtoPack(data, UniqueParObjectId());

  // add base class Element
  DRT::Element::Pack(data);

  AddtoPack(data, (int)celltype_);

  AddToPack(data, solid_ele_property_);

  AddtoPack(data, poro_ele_property_.porotype);

  AddtoPack(data, poro_ele_property_.impltype);

  data.AddtoPack(material_post_setup_);

  // anisotropic_permeability_directions_
  auto size = static_cast<int>(anisotropic_permeability_directions_.size());
  AddtoPack(data, size);  // TODO:: necessary?
  for (int i = 0; i < size; ++i) AddtoPack(data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = static_cast<int>(anisotropic_permeability_nodal_coeffs_.size());
  AddtoPack(data, size);  // TODO:: necessary?
  for (int i = 0; i < size; ++i) AddtoPack(data, anisotropic_permeability_nodal_coeffs_[i]);

  // optional data, e.g., EAS data
  DRT::ELEMENTS::Pack(solid_calc_variant_, data);
  DRT::ELEMENTS::Pack(solidporo_calc_variant_, data);
}

void DRT::ELEMENTS::SolidPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  if (ExtractInt(position, data) != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Element::Unpack(basedata);

  celltype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));

  ExtractFromPack(position, data, solid_ele_property_);

  poro_ele_property_.porotype = static_cast<INPAR::PORO::PoroType>(ExtractInt(position, data));

  poro_ele_property_.impltype = static_cast<INPAR::SCATRA::ImplType>(ExtractInt(position, data));

  CORE::COMM::ParObject::ExtractfromPack(position, data, material_post_setup_);

  // anisotropic_permeability_directions_
  int size = 0;
  ExtractfromPack(position, data, size);
  anisotropic_permeability_directions_.resize(size, std::vector<double>(3, 0.0));
  for (int i = 0; i < size; ++i)
    ExtractfromPack(position, data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = 0;
  ExtractfromPack(position, data, size);
  anisotropic_permeability_nodal_coeffs_.resize(size, std::vector<double>(this->NumNode(), 0.0));
  for (int i = 0; i < size; ++i)
    ExtractfromPack(position, data, anisotropic_permeability_nodal_coeffs_[i]);

  // reset solid and poro interfaces
  solid_calc_variant_ = CreateSolidCalculationInterface(celltype_, solid_ele_property_);
  solidporo_calc_variant_ = CreateSolidPoroCalculationInterface(*this, GetElePoroType());

  DRT::ELEMENTS::Unpack(solid_calc_variant_, position, data);
  DRT::ELEMENTS::Unpack(solidporo_calc_variant_, position, data);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

void DRT::ELEMENTS::SolidPoro::VisNames(std::map<std::string, int>& names)
{
  DRT::Element::VisNames(names);
  SolidPoroMaterial().VisNames(names);
}

bool DRT::ELEMENTS::SolidPoro::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return SolidPoroMaterial().VisData(name, data, Id());
}

MAT::StructPoro& DRT::ELEMENTS::SolidPoro::StructPoroMaterial(int nummat) const
{
  auto porostruct_mat =
      Teuchos::rcp_dynamic_cast<MAT::StructPoro>(DRT::Element::Material(nummat), true);

  if (porostruct_mat == Teuchos::null) dserror("cast to poro material failed");

  if (porostruct_mat->MaterialType() != INPAR::MAT::m_structporo and
      porostruct_mat->MaterialType() != INPAR::MAT::m_structpororeaction and
      porostruct_mat->MaterialType() != INPAR::MAT::m_structpororeactionECM)
    dserror("invalid structure material for poroelasticity");

  return *porostruct_mat;
}
MAT::FluidPoroMultiPhase& DRT::ELEMENTS::SolidPoro::FluidPoroMultiMaterial(int nummat) const
{
  if (this->NumMaterial() <= 1)
  {
    dserror("No second material defined for SolidPoro element %i", Id());
  }

  auto fluidmulti_mat =
      Teuchos::rcp_dynamic_cast<MAT::FluidPoroMultiPhase>(DRT::Element::Material(1), true);

  if (fluidmulti_mat == Teuchos::null) dserror("cast to multiphase fluid poro material failed");
  if (fluidmulti_mat->MaterialType() != INPAR::MAT::m_fluidporo_multiphase and
      fluidmulti_mat->MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("invalid fluid material for poro-multiphase-elasticity");
  if (fluidmulti_mat->NumFluidPhases() == 0)
  {
    dserror(
        "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE = 0 currently not supported since this requires "
        "an adaption of the definition of the solid pressure");
  }
  return *fluidmulti_mat;
}

BACI_NAMESPACE_CLOSE
