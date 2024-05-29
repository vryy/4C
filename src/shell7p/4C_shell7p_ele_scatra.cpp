/*! \file
\brief A 2D shell element with ScaTra functionality

\level 3
*/

#include "4C_shell7p_ele_scatra.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_shell7p_ele_factory.hpp"
#include "4C_shell7p_ele_interface_serializable.hpp"
#include "4C_shell7p_line.hpp"
#include "4C_shell7p_utils.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <typename Interface>
  void TryPackInterface(const Interface& interface, CORE::COMM::PackBuffer& data)
  {
    std::shared_ptr<DRT::ELEMENTS::SHELL::Serializable> serializable_interface =
        std::dynamic_pointer_cast<DRT::ELEMENTS::SHELL::Serializable>(interface);
    if (serializable_interface != nullptr) serializable_interface->Pack(data);
  }

  template <typename Interface>
  void TryUnpackInterface(
      Interface& interface, std::size_t& position, const std::vector<char>& data)
  {
    std::shared_ptr<DRT::ELEMENTS::SHELL::Serializable> serializable_shell_interface =
        std::dynamic_pointer_cast<DRT::ELEMENTS::SHELL::Serializable>(interface);
    if (serializable_shell_interface != nullptr)
      serializable_shell_interface->Unpack(position, data);
  }

}  // namespace

DRT::ELEMENTS::Shell7pScatraType DRT::ELEMENTS::Shell7pScatraType::instance_;


DRT::ELEMENTS::Shell7pScatraType& DRT::ELEMENTS::Shell7pScatraType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::Shell7pScatraType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Shell7pScatra(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::Shell7pScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SHELL7PSCATRA") return Create(id, owner);
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::Shell7pScatraType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::Shell7pScatra(id, owner));
}

void DRT::ELEMENTS::Shell7pScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defsgeneral = definitions["SHELL7PSCATRA"];

  defsgeneral["QUAD4"] = INPUT::LineDefinition::Builder()
                             .AddIntVector("QUAD4", 4)
                             .AddNamedInt("MAT")
                             .AddNamedDouble("THICK")
                             .AddNamedString("EAS")
                             .AddString("EAS2")
                             .AddString("EAS3")
                             .AddString("EAS4")
                             .AddString("EAS5")
                             .AddNamedDouble("SDC")
                             .AddOptionalTag("ANS")
                             .add_optional_named_double_vector("RAD", 3)
                             .add_optional_named_double_vector("AXI", 3)
                             .add_optional_named_double_vector("CIR", 3)
                             .add_optional_named_double_vector("FIBER1", 3)
                             .add_optional_named_double_vector("FIBER2", 3)
                             .add_optional_named_double_vector("FIBER3", 3)
                             .add_optional_named_string("TYPE")
                             .Build();

  defsgeneral["QUAD8"] = INPUT::LineDefinition::Builder()
                             .AddIntVector("QUAD8", 8)
                             .AddNamedInt("MAT")
                             .AddNamedDouble("THICK")
                             .AddNamedString("EAS")
                             .AddString("EAS2")
                             .AddString("EAS3")
                             .AddString("EAS4")
                             .AddString("EAS5")
                             .AddNamedDouble("SDC")
                             .AddOptionalTag("ANS")
                             .add_optional_named_double_vector("RAD", 3)
                             .add_optional_named_double_vector("AXI", 3)
                             .add_optional_named_double_vector("CIR", 3)
                             .add_optional_named_double_vector("FIBER1", 3)
                             .add_optional_named_double_vector("FIBER2", 3)
                             .add_optional_named_double_vector("FIBER3", 3)
                             .add_optional_named_string("TYPE")
                             .Build();

  defsgeneral["QUAD9"] = INPUT::LineDefinition::Builder()
                             .AddIntVector("QUAD9", 9)
                             .AddNamedInt("MAT")
                             .AddNamedDouble("THICK")
                             .AddNamedString("EAS")
                             .AddString("EAS2")
                             .AddString("EAS3")
                             .AddString("EAS4")
                             .AddString("EAS5")
                             .AddNamedDouble("SDC")
                             .AddOptionalTag("ANS")
                             .add_optional_named_double_vector("RAD", 3)
                             .add_optional_named_double_vector("AXI", 3)
                             .add_optional_named_double_vector("CIR", 3)
                             .add_optional_named_double_vector("FIBER1", 3)
                             .add_optional_named_double_vector("FIBER2", 3)
                             .add_optional_named_double_vector("FIBER3", 3)
                             .add_optional_named_string("TYPE")
                             .Build();

  defsgeneral["TRI3"] = INPUT::LineDefinition::Builder()
                            .AddIntVector("TRI3", 3)
                            .AddNamedInt("MAT")
                            .AddNamedDouble("THICK")
                            .AddNamedDouble("SDC")
                            .add_optional_named_double_vector("RAD", 3)
                            .add_optional_named_double_vector("AXI", 3)
                            .add_optional_named_double_vector("CIR", 3)
                            .add_optional_named_double_vector("FIBER1", 3)
                            .add_optional_named_double_vector("FIBER2", 3)
                            .add_optional_named_double_vector("FIBER3", 3)
                            .add_optional_named_string("TYPE")
                            .Build();

  defsgeneral["TRI6"] = INPUT::LineDefinition::Builder()
                            .AddIntVector("TRI6", 6)
                            .AddNamedInt("MAT")
                            .AddNamedDouble("THICK")
                            .AddNamedDouble("SDC")
                            .add_optional_named_double_vector("RAD", 3)
                            .add_optional_named_double_vector("AXI", 3)
                            .add_optional_named_double_vector("CIR", 3)
                            .add_optional_named_double_vector("FIBER1", 3)
                            .add_optional_named_double_vector("FIBER2", 3)
                            .add_optional_named_double_vector("FIBER3", 3)
                            .add_optional_named_string("TYPE")
                            .Build();
}

int DRT::ELEMENTS::Shell7pScatraType::Initialize(DRT::Discretization& dis)
{
  STR::UTILS::SHELL::DIRECTOR::SetupShellElementDirectors(*this, dis);

  return 0;
}



CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::Shell7pScatraType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  auto* shell = dynamic_cast<DRT::ELEMENTS::Shell7pScatra*>(node.Elements()[0]);
  if (!shell) FOUR_C_THROW("Cannot cast to Shell");
  int j;
  for (j = 0; j < shell->num_node(); ++j)
    if (shell->Nodes()[j]->Id() == node.Id()) break;
  if (j == shell->num_node()) FOUR_C_THROW("Can't find matching node..!");
  double half_thickness = shell->GetThickness() / 2.0;

  // set director
  const CORE::LINALG::SerialDenseMatrix nodal_directors = shell->GetDirectors();
  CORE::LINALG::Matrix<SHELL::DETAIL::num_dim, 1> director(true);
  for (int dim = 0; dim < SHELL::DETAIL::num_dim; ++dim)
    director(dim, 0) = nodal_directors(j, dim) * half_thickness;

  return STR::UTILS::SHELL::ComputeShellNullSpace(node, x0, director);
}

void DRT::ELEMENTS::Shell7pScatraType::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  STR::UTILS::SHELL::NodalBlockInformationShell(dwele, numdf, dimns, nv, np);
}


DRT::ELEMENTS::Shell7pScatra::Shell7pScatra(const DRT::ELEMENTS::Shell7pScatra& other)
    : CORE::Elements::Element(other),
      distype_(other.distype_),
      interface_ptr_(other.interface_ptr_),
      eletech_(other.eletech_),
      thickness_(other.thickness_),
      nodal_directors_(other.nodal_directors_),
      material_post_setup_(other.material_post_setup_),
      impltype_(other.impltype_)
{
  // reset shell calculation interface
  shell_interface_ = Shell7pFactory::provide_shell7p_calculation_interface(other, other.eletech_);
}


DRT::ELEMENTS::Shell7pScatra& DRT::ELEMENTS::Shell7pScatra::operator=(
    const DRT::ELEMENTS::Shell7pScatra& other)
{
  if (this == &other) return *this;
  CORE::Elements::Element::operator=(other);
  distype_ = other.distype_;
  interface_ptr_ = other.interface_ptr_;
  eletech_ = other.eletech_;
  thickness_ = other.thickness_;
  nodal_directors_ = other.nodal_directors_;
  material_post_setup_ = other.material_post_setup_;
  impltype_ = other.impltype_;

  shell_interface_ = Shell7pFactory::provide_shell7p_calculation_interface(other, other.eletech_);
  return *this;
}

CORE::Elements::Element* DRT::ELEMENTS::Shell7pScatra::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::Shell7pScatra(*this);
  return newelement;
}

void DRT::ELEMENTS::Shell7pScatra::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  CORE::Elements::Element::Pack(data);
  // discretization type
  AddtoPack(data, (int)distype_);
  // element technology
  AddtoPack(data, eletech_);
  // thickness in reference frame
  AddtoPack(data, thickness_);
  // nodal_directors
  AddtoPack(data, nodal_directors_);
  // Setup flag for material post setup
  data.AddtoPack(material_post_setup_);
  // pack impltype
  AddtoPack(data, impltype_);
  // optional data, e.g., EAS data, current thickness,..
  TryPackInterface(shell_interface_, data);
}


void DRT::ELEMENTS::Shell7pScatra::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // discretization type
  distype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));
  // element technology
  ExtractfromPack(position, data, eletech_);
  // thickness in reference frame
  ExtractfromPack(position, data, thickness_);
  // nodal director
  ExtractfromPack(position, data, nodal_directors_);
  // Setup flag for material post setup
  CORE::COMM::ParObject::ExtractfromPack(position, data, material_post_setup_);
  // extract impltype
  impltype_ =
      static_cast<INPAR::SCATRA::ImplType>(CORE::COMM::ParObject::ExtractInt(position, data));
  // reset shell calculation interface
  shell_interface_ = Shell7pFactory::provide_shell7p_calculation_interface(*this, eletech_);

  TryUnpackInterface(shell_interface_, position, data);
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

Teuchos::RCP<MAT::So3Material> DRT::ELEMENTS::Shell7pScatra::SolidMaterial(int nummat) const
{
  return Teuchos::rcp_dynamic_cast<MAT::So3Material>(
      CORE::Elements::Element::Material(nummat), true);
}

void DRT::ELEMENTS::Shell7pScatra::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<CORE::Elements::ParamsInterface>>("interface"));
  }
  else
  {
    interface_ptr_ = Teuchos::null;
  }
}


void DRT::ELEMENTS::Shell7pScatra::VisNames(std::map<std::string, int>& names)
{
  std::string result_thickness = "thickness";
  names[result_thickness] = 1;
  SolidMaterial()->VisNames(names);
}  // VisNames()


bool DRT::ELEMENTS::Shell7pScatra::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (CORE::Elements::Element::VisData(name, data)) return true;

  shell_interface_->VisData(name, data);

  return SolidMaterial()->VisData(name, data, Id());

}  // VisData()


void DRT::ELEMENTS::Shell7pScatra::Print(std::ostream& os) const
{
  os << "Shell7pScatra ";
  os << " discretization type: " << CORE::FE::CellTypeToString(distype_).c_str();
  Element::Print(os);
}

std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::Shell7pScatra::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<Shell7pLine, Shell7pScatra>(
      CORE::COMM::buildLines, *this);
}

std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::Shell7pScatra::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

int DRT::ELEMENTS::Shell7pScatra::NumLine() const
{
  return CORE::FE::getNumberOfElementLines(distype_);
}


int DRT::ELEMENTS::Shell7pScatra::NumSurface() const { return 1; }


bool DRT::ELEMENTS::Shell7pScatra::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  STR::ELEMENTS::ShellData shell_data = {};

  // set discretization type
  distype_ = CORE::FE::StringToCellType(distype);

  // set thickness in reference frame
  linedef->ExtractDouble("THICK", thickness_);
  if (thickness_ <= 0) FOUR_C_THROW("Shell element thickness needs to be > 0");
  shell_data.thickness = thickness_;

  // extract number of EAS parameters for different locking types
  STR::ELEMENTS::ShellLockingTypes locking_types = {};
  if (linedef->HaveNamed("EAS"))
  {
    eletech_.insert(INPAR::STR::EleTech::eas);
    STR::UTILS::SHELL::READELEMENT::ReadAndSetLockingTypes(distype_, linedef, locking_types);
  }

  // set calculation interface pointer
  shell_interface_ = Shell7pFactory::provide_shell7p_calculation_interface(*this, eletech_);

  // read and set ANS technology for element
  if (linedef->HaveNamed("ANS"))
  {
    shell_data.num_ans = STR::UTILS::SHELL::READELEMENT::ReadAndSetNumANS(distype_);
  }
  // read SDC
  linedef->ExtractDouble("SDC", shell_data.sdc);

  // read and set number of material model
  SetMaterial(0, MAT::Factory(STR::UTILS::SHELL::READELEMENT::ReadAndSetElementMaterial(linedef)));

  // setup shell calculation interface
  shell_interface_->Setup(*this, *SolidMaterial(), linedef, locking_types, shell_data);
  if (!material_post_setup_)
  {
    shell_interface_->material_post_setup(*this, *SolidMaterial());
    material_post_setup_ = true;
  }
  // read implementation type for scatra
  std::string impltype;
  linedef->ExtractString("TYPE", impltype);

  if (impltype == "Undefined")
    impltype_ = INPAR::SCATRA::impltype_undefined;
  else if (impltype == "AdvReac")
    impltype_ = INPAR::SCATRA::impltype_advreac;
  else if (impltype == "CardMono")
    impltype_ = INPAR::SCATRA::impltype_cardiac_monodomain;
  else if (impltype == "Chemo")
    impltype_ = INPAR::SCATRA::impltype_chemo;
  else if (impltype == "ChemoReac")
    impltype_ = INPAR::SCATRA::impltype_chemoreac;
  else if (impltype == "Loma")
    impltype_ = INPAR::SCATRA::impltype_loma;
  else if (impltype == "RefConcReac")
    impltype_ = INPAR::SCATRA::impltype_refconcreac;
  else if (impltype == "Std")
    impltype_ = INPAR::SCATRA::impltype_std;
  else
    FOUR_C_THROW("Invalid implementation type for Shell7pScatra elements!");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
