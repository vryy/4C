/*! \file
\brief

\level 1

*----------------------------------------------------------------------*/

#include "solid_ele.H"
#include <memory>
#include "solid_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_so3/so_line.H"
#include "../drt_so3/so_surface.H"
#include "../linalg/linalg_utils_nullspace.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"
#include "../drt_mat/so3_material.H"
#include "solid_ele_factory.H"
#include "solid_ele_interface.H"


DRT::ELEMENTS::SolidType DRT::ELEMENTS::SolidType::instance_;

DRT::ELEMENTS::SolidType& DRT::ELEMENTS::SolidType::Instance() { return instance_; }

void DRT::ELEMENTS::SolidType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defsgeneral = definitions["SOLID"];

  defsgeneral["HEX8"]
      .AddIntVector("HEX8", 8)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedString("EAS")
      .AddOptionalNamedString("FBAR")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["HEX18"]
      .AddIntVector("HEX18", 18)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["HEX20"]
      .AddIntVector("HEX20", 20)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["HEX27"]
      .AddIntVector("HEX27", 27)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["TET4"]
      .AddIntVector("TET4", 4)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["TET10"]
      .AddIntVector("TET10", 10)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["WEDGE6"]
      .AddIntVector("WEDGE6", 6)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["WEDGE15"]
      .AddIntVector("WEDGE15", 15)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["PYRAMID5"]
      .AddIntVector("PYRAMID5", 5)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["NURBS8"]
      .AddIntVector("NURBS8", 8)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);

  defsgeneral["NURBS27"]
      .AddIntVector("NURBS27", 27)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3);
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
  // todo: to this combined for 2D and 3D
  numdf = 3;
  dimns = 6;

  nv = 3;
}

Teuchos::SerialDenseMatrix<int, double> DRT::ELEMENTS::SolidType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  // todo make this work for 3D
  dserror("Not implemented yet!");
  exit(1);
}

DRT::ELEMENTS::Solid::Solid(int id, int owner)
    : DRT::Element(id, owner),
      distype_(DRT::Element::dis_none),
      kintype_(INPAR::STR::kinem_vague),
      eastype_(STR::ELEMENTS::EASType::eastype_undefined),
      interface_ptr_(Teuchos::null)
{
}

DRT::ELEMENTS::Solid::Solid(const DRT::ELEMENTS::Solid& old)
    : DRT::Element(old),
      distype_(old.distype_),
      kintype_(old.kintype_),
      eastype_(old.eastype_),
      interface_ptr_(old.interface_ptr_)
{
}

DRT::Element* DRT::ELEMENTS::Solid::Clone() const
{
  return nullptr;
  //  return new DRT::ELEMENTS::Solid(*this);
}

int DRT::ELEMENTS::Solid::NumLine() const { return DRT::UTILS::getNumberOfElementLines(distype_); }

int DRT::ELEMENTS::Solid::NumSurface() const
{
  return DRT::UTILS::getNumberOfElementSurfaces(distype_);
}

int DRT::ELEMENTS::Solid::NumVolume() const
{
  return DRT::UTILS::getNumberOfElementVolumes(distype_);
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Solid::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumLine() > 1)  // 1D boundary element and 2D/3D parent element
  {
    return DRT::UTILS::ElementBoundaryFactory<StructuralLine, Solid>(DRT::UTILS::buildLines, this);
  }
  else if (NumLine() ==
           1)  // 1D boundary element and 1D parent element -> body load (calculated in evaluate)
  {
    // 1D (we return the element itself)
    std::vector<Teuchos::RCP<Element>> surfaces(1);
    surfaces[0] = Teuchos::rcp(this, false);
    return surfaces;
  }
  else
  {
    dserror("Lines() does not exist for points ");
    return DRT::Element::Surfaces();
  }
}

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Solid::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:

  if (NumSurface() > 1)  // 2D boundary element and 3D parent element
    return DRT::UTILS::ElementBoundaryFactory<StructuralSurface, Solid>(
        DRT::UTILS::buildSurfaces, this);
  else if (NumSurface() ==
           1)  // 2D boundary element and 2D parent element -> body load (calculated in evaluate)
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

std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Solid::Volumes()
{
  if (NumVolume() ==
      1)  // 3D boundary element and a 3D parent element -> body load (calculated in evaluate)
  {
    std::vector<Teuchos::RCP<Element>> volumes(1);
    volumes[0] = Teuchos::rcp(this, false);
    return volumes;
  }
  else  //
  {
    dserror("Volumes() does not exist for 1D/2D-elements");
    return DRT::Element::Surfaces();
  }
}

void DRT::ELEMENTS::Solid::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // Object Id
  AddtoPack(data, UniqueParObjectId());

  // add base class Element
  DRT::Element::Pack(data);

  // discretization type
  AddtoPack(data, (int)distype_);

  // kinematic type
  AddtoPack(data, (int)kintype_);

  // element technology7
  AddtoPack(data, eletech_);

  // eas type
  AddtoPack(data, eastype_);

  // Setup flag
  data.AddtoPack(material_post_setup_);
}

void DRT::ELEMENTS::Solid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  if (ExtractInt(position, data) != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Element::Unpack(basedata);

  // discretization type
  distype_ = static_cast<DRT::Element::DiscretizationType>(ExtractInt(position, data));
  // kinematic type
  kintype_ = static_cast<INPAR::STR::KinemType>(ExtractInt(position, data));
  // element technology
  DRT::ParObject::ExtractfromPack(position, data, eletech_);
  // eas type
  eastype_ = static_cast<::STR::ELEMENTS::EASType>(ExtractInt(position, data));
  // Setup flag
  DRT::ParObject::ExtractfromPack(position, data, material_post_setup_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

void DRT::ELEMENTS::Solid::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  // if (interface_ptr_ != Teuchos::null) return;
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
  distype_ = DRT::StringToDistype(distype);
  solid_interface_ = SolidFactory::ProvideImpl(this);

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);
  //  SolidMaterial()->Setup(STR::UTILS::DisTypeToNgpOptGaussRule(Shape()), linedef);

  // kinematic type
  std::string kinem;
  linedef->ExtractString("KINEM", kinem);
  if (kinem == "nonlinear")
    SetKinematicType(INPAR::STR::kinem_nonlinearTotLag);
  else if (kinem == "linear")
    SetKinematicType(INPAR::STR::kinem_linear);
  else
    dserror("unknown kinematic type %s", kinem.c_str());

  if (linedef->HaveNamed("EAS"))
  {
    std::string eastype;
    linedef->ExtractString("EAS", eastype);
    if (Shape() == DRT::Element::hex8)
    {
      if (eastype == "mild")
      {
        eastype_ = ::STR::ELEMENTS::EASType::eastype_h8_9;
        eletech_.insert(INPAR::STR::EleTech::eas);
      }
      else if (eastype == "full")
      {
        eastype_ = ::STR::ELEMENTS::EASType::eastype_h8_21;
        eletech_.insert(INPAR::STR::EleTech::eas);
      }
      else if (eastype == "none")
      {
        eastype_ = ::STR::ELEMENTS::EASType::soh8_easnone;
      }
      else
        dserror("unrecognized eas type for hex8: %s", eastype.c_str());
    }
    else
      dserror("no EAS allowed for this element shape");
  }

  DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->Setup(this, linedef);
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
}  // VisNames()

bool DRT::ELEMENTS::Solid::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, Id());

}  // VisData()
int DRT::ELEMENTS::Solid::PostEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  return 0;
}
