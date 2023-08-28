/*! \file

\brief Implementation of the solid-poro element

\level 1
*/

#include "baci_solid_ele_poro.H"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_mat_fluidporo_multiphase.H"
#include "baci_mat_structporo.H"
#include "baci_so3_line.H"
#include "baci_so3_nullspace.H"
#include "baci_so3_surface.H"
#include "baci_solid_ele_factory.H"
#include "baci_solid_ele_interface_serializable.H"
#include "baci_solid_ele_poro_factory.H"
#include "baci_solid_ele_poro_utils.H"
#include "baci_solid_ele_utils.H"

#include <memory>

namespace
{
  template <typename Interface>
  void TryPackInterface(const Interface& interface, DRT::PackBuffer& data)
  {
    std::shared_ptr<DRT::ELEMENTS::Serializable> serializable_interface =
        std::dynamic_pointer_cast<DRT::ELEMENTS::Serializable>(interface);
    if (serializable_interface != nullptr) serializable_interface->Pack(data);
  }

  template <typename Interface>
  void TryUnpackInterface(
      Interface& interface, std::size_t& position, const std::vector<char>& data)
  {
    std::shared_ptr<DRT::ELEMENTS::Serializable> serializable_solid_interface =
        std::dynamic_pointer_cast<DRT::ELEMENTS::Serializable>(interface);
    if (serializable_solid_interface != nullptr)
      serializable_solid_interface->Unpack(position, data);
  }
}  // namespace


DRT::ELEMENTS::SolidPoroType DRT::ELEMENTS::SolidPoroType::instance_;

DRT::ELEMENTS::SolidPoroType& DRT::ELEMENTS::SolidPoroType::Instance() { return instance_; }

void DRT::ELEMENTS::SolidPoroType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defsgeneral = definitions["SOLIDPORO"];

  defsgeneral["HEX8"] = INPUT::LineDefinition::Builder()
                            .AddIntVector("HEX8", 8)
                            .AddNamedInt("MAT")
                            .AddNamedString("KINEM")
                            .AddOptionalNamedString("EAS")
                            .AddOptionalTag("FBAR")
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
                            .AddOptionalNamedDoubleVector("POROANISODIR3", 3)
                            .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS1", 8)
                            .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS2", 8)
                            .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS3", 8)
                            .Build();

  defsgeneral["HEX27"] = INPUT::LineDefinition::Builder()
                             .AddIntVector("HEX27", 27)
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
                             .AddOptionalNamedDoubleVector("POROANISODIR3", 3)
                             .Build();

  defsgeneral["TET4"] = INPUT::LineDefinition::Builder()
                            .AddIntVector("TET4", 4)
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
                            .AddOptionalNamedDoubleVector("POROANISODIR3", 3)
                            .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS1", 8)
                            .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS2", 8)
                            .AddOptionalNamedDoubleVector("POROANISONODALCOEFFS3", 8)
                            .Build();

  defsgeneral["TET10"] = INPUT::LineDefinition::Builder()
                             .AddIntVector("TET10", 10)
                             .AddNamedInt("MAT")
                             .AddNamedString("KINEM")
                             .AddOptionalNamedDoubleVector("RAD", 3)
                             .AddOptionalNamedDoubleVector("AXI", 3)
                             .AddOptionalNamedDoubleVector("CIR", 3)
                             .AddOptionalNamedDoubleVector("FIBER1", 3)
                             .AddOptionalNamedDoubleVector("FIBER2", 3)
                             .AddOptionalNamedDoubleVector("FIBER3", 3)
                             .AddOptionalNamedString("POROTYPE")
                             .AddOptionalNamedString("TYPE")
                             .AddOptionalNamedDoubleVector("POROANISODIR1", 3)
                             .AddOptionalNamedDoubleVector("POROANISODIR2", 3)
                             .AddOptionalNamedDoubleVector("POROANISODIR3", 3)
                             .Build();
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidPoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDPORO") return Create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SolidPoroType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::SolidPoro(id, owner));
}

DRT::ParObject* DRT::ELEMENTS::SolidPoroType::Create(const std::vector<char>& data)
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

DRT::ELEMENTS::SolidPoro::SolidPoro(int id, int owner)
    : DRT::Element(id, owner),
      distype_(DRT::Element::dis_none),
      kintype_(INPAR::STR::kinem_vague),
      eastype_(STR::ELEMENTS::EasType::eastype_undefined),
      porotype_(INPAR::PORO::PoroType::undefined),
      impltype_(INPAR::SCATRA::impltype_undefined)
{
}

DRT::ELEMENTS::SolidPoro::SolidPoro(const DRT::ELEMENTS::SolidPoro& other)
    : DRT::Element(other),
      distype_(other.distype_),
      kintype_(other.kintype_),
      eastype_(other.eastype_),
      porotype_(other.porotype_),
      impltype_(other.impltype_),
      anisotropic_permeability_directions_(other.anisotropic_permeability_directions_),
      anisotropic_permeability_nodal_coeffs_(other.anisotropic_permeability_nodal_coeffs_),
      interface_ptr_(other.interface_ptr_),
      material_post_setup_(other.material_post_setup_)
{
  // create own solid and poro interface on copy
  CreateSolidCalculationInterface(*this, GetEleTech(), GetEleKinematicType(), GetEAStype());
  solidporo_interface_ = CreateSolidPoroCalculationInterface(*this, GetElePoroType());
}

DRT::ELEMENTS::SolidPoro& ::DRT::ELEMENTS::SolidPoro::operator=(
    const DRT::ELEMENTS::SolidPoro& other)
{
  distype_ = other.distype_;
  kintype_ = other.kintype_;
  eastype_ = other.eastype_;
  porotype_ = other.porotype_;
  impltype_ = other.impltype_;
  anisotropic_permeability_directions_ = other.anisotropic_permeability_directions_;
  anisotropic_permeability_nodal_coeffs_ = other.anisotropic_permeability_nodal_coeffs_;
  interface_ptr_ = other.interface_ptr_;
  material_post_setup_ = other.material_post_setup_;

  // create own solid and poro interface on copy assignment
  solid_interface_ =
      CreateSolidCalculationInterface(*this, GetEleTech(), GetEleKinematicType(), GetEAStype());
  solidporo_interface_ = CreateSolidPoroCalculationInterface(*this, GetElePoroType());

  return *this;
}

DRT::Element* DRT::ELEMENTS::SolidPoro::Clone() const
{
  return new DRT::ELEMENTS::SolidPoro(*this);
}

int DRT::ELEMENTS::SolidPoro::NumLine() const
{
  return CORE::DRT::UTILS::getNumberOfElementLines(distype_);
}

int DRT::ELEMENTS::SolidPoro::NumSurface() const
{
  return CORE::DRT::UTILS::getNumberOfElementSurfaces(distype_);
}

int DRT::ELEMENTS::SolidPoro::NumVolume() const
{
  return CORE::DRT::UTILS::getNumberOfElementVolumes(distype_);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::SolidPoro::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine() > 1)  // 1D boundary element and 2D/3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<StructuralLine, SolidPoro>(
        DRT::UTILS::buildLines, this);
  }
  else if (NumLine() == 1)  // 1D boundary element and 1D parent element -> body load
                            // (calculated in evaluate)
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element>> lines(1);
    lines[0] = Teuchos::rcp(this, false);
    return lines;
  }
  else
  {
    dserror("Lines() does not exist for points ");
    return DRT::Element::Lines();
  }
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::SolidPoro::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumSurface() > 1)  // 2D boundary element and 3D parent element
    return DRT::UTILS::ElementBoundaryFactory<StructuralSurface, SolidPoro>(
        DRT::UTILS::buildSurfaces, this);
  else if (NumSurface() == 1)  // 2D boundary element and 2D parent element -> body load
                               // (calculated in evaluate)
  {
    // 2D (we return the element itself)
    std::vector<Teuchos::RCP<Element>> surfaces(1);
    surfaces[0] = Teuchos::rcp(this, false);
    return surfaces;
  }
  else  // 1D elements
  {
    dserror("Surfaces() does not exist for 1D-element ");
    return DRT::Element::Surfaces();
  }
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::SolidPoro::Volumes()
{
  if (NumVolume() ==
      1)  // 3D boundary element and a 3D parent element -> body load (calculated in evaluate)
  {
    std::vector<Teuchos::RCP<Element>> volumes(1);
    volumes[0] = Teuchos::rcp(this, false);
    return volumes;
  }
  else
  {
    dserror("Volumes() does not exist for 1D/2D-elements");
    return DRT::Element::Volumes();
  }
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
    const std::string& eletype, const std::string& eledistype, DRT::INPUT::LineDefinition* linedef)
{
  // read base element
  // set discretization type
  distype_ = DRT::StringToDistype(eledistype);
  anisotropic_permeability_directions_.resize(3, std::vector<double>(3, 0.0));
  anisotropic_permeability_nodal_coeffs_.resize(3, std::vector<double>(this->NumNode(), 0.0));

  // read number of material model
  SetMaterial(STR::UTILS::READELEMENT::ReadElementMaterial(linedef));

  // kinematic type
  SetKinematicType(STR::UTILS::READELEMENT::ReadElementKinematicType(linedef));


  if (linedef->HaveNamed("EAS"))
  {
    if (Shape() == DRT::Element::hex8)
    {
      STR::UTILS::READELEMENT::ReadAndSetEAS(linedef, eastype_, eletech_);
    }
    else
      dserror("no EAS allowed for this element shape");
  }

  // read scalar transport implementation type
  if (linedef->HaveNamed("POROTYPE"))
  {
    porotype_ = STR::UTILS::READELEMENT::ReadPoroType(linedef);
  }
  else
  {
    porotype_ = INPAR::PORO::PoroType::undefined;
  }

  // read scalar transport implementation type
  if (linedef->HaveNamed("TYPE"))
  {
    impltype_ = STR::UTILS::READELEMENT::ReadType(linedef);
  }
  else
  {
    impltype_ = INPAR::SCATRA::impltype_undefined;
  }

  ReadAnisotropicPermeabilityDirectionsFromElementLineDefinition(linedef);
  ReadAnisotropicPermeabilityNodalCoeffsFromElementLineDefinition(linedef);

  solid_interface_ =
      CreateSolidCalculationInterface(*this, GetEleTech(), GetEleKinematicType(), GetEAStype());
  solidporo_interface_ = CreateSolidPoroCalculationInterface(*this, GetElePoroType());

  // setup solid material
  solid_interface_->Setup(StructPoroMaterial(), linedef);

  // setup poro material
  solidporo_interface_->PoroSetup(StructPoroMaterial(), linedef);

  return true;
}

void DRT::ELEMENTS::SolidPoro::ReadAnisotropicPermeabilityDirectionsFromElementLineDefinition(
    DRT::INPUT::LineDefinition* linedef)
{
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string definition_name = "POROANISODIR" + std::to_string(dim + 1);
    if (linedef->HaveNamed(definition_name))
      linedef->ExtractDoubleVector(definition_name, anisotropic_permeability_directions_[dim]);
  }
}

void DRT::ELEMENTS::SolidPoro::ReadAnisotropicPermeabilityNodalCoeffsFromElementLineDefinition(
    DRT::INPUT::LineDefinition* linedef)
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

void DRT::ELEMENTS::SolidPoro::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  AddtoPack(data, UniqueParObjectId());

  // add base class Element
  DRT::Element::Pack(data);

  AddtoPack(data, (int)distype_);

  AddtoPack(data, (int)kintype_);

  AddtoPack(data, eletech_);

  AddtoPack(data, eastype_);

  AddtoPack(data, porotype_);

  AddtoPack(data, impltype_);

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
  TryPackInterface(solid_interface_, data);
  TryPackInterface(solidporo_interface_, data);
}

void DRT::ELEMENTS::SolidPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  if (ExtractInt(position, data) != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Element::Unpack(basedata);

  distype_ = static_cast<DRT::Element::DiscretizationType>(ExtractInt(position, data));

  kintype_ = static_cast<INPAR::STR::KinemType>(ExtractInt(position, data));

  DRT::ParObject::ExtractfromPack(position, data, eletech_);

  eastype_ = static_cast<::STR::ELEMENTS::EasType>(ExtractInt(position, data));

  porotype_ = static_cast<INPAR::PORO::PoroType>(ExtractInt(position, data));

  impltype_ = static_cast<INPAR::SCATRA::ImplType>(ExtractInt(position, data));

  DRT::ParObject::ExtractfromPack(position, data, material_post_setup_);

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
  solid_interface_ =
      CreateSolidCalculationInterface(*this, GetEleTech(), GetEleKinematicType(), GetEAStype());
  solidporo_interface_ = CreateSolidPoroCalculationInterface(*this, GetElePoroType());

  TryUnpackInterface(solid_interface_, position, data);
  TryUnpackInterface(solidporo_interface_, position, data);

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
