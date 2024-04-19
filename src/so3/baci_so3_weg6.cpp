/*----------------------------------------------------------------------*/
/*! \file

\brief Solid Wedge6 Element

\level 1


*----------------------------------------------------------------------*/

#include "baci_so3_weg6.hpp"

#include "baci_comm_utils_factory.hpp"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "baci_global_data.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_lib_discret.hpp"
#include "baci_mat_so3_material.hpp"
#include "baci_so3_line.hpp"
#include "baci_so3_nullspace.hpp"
#include "baci_so3_prestress.hpp"
#include "baci_so3_prestress_service.hpp"
#include "baci_so3_surface.hpp"
#include "baci_so3_utils.hpp"
#include "baci_utils_exceptions.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::SoWeg6Type DRT::ELEMENTS::SoWeg6Type::instance_;

DRT::ELEMENTS::SoWeg6Type& DRT::ELEMENTS::SoWeg6Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoWeg6Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SoWeg6(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoWeg6Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoWeg6(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoWeg6Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoWeg6(id, owner));
  return ele;
}


void DRT::ELEMENTS::SoWeg6Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SoWeg6Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::SoWeg6Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["WEDGE6"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("WEDGE6", 6)
                       .AddNamedInt("MAT")
                       .AddNamedString("KINEM")
                       .AddOptionalNamedDoubleVector("RAD", 3)
                       .AddOptionalNamedDoubleVector("AXI", 3)
                       .AddOptionalNamedDoubleVector("CIR", 3)
                       .AddOptionalNamedDoubleVector("FIBER1", 3)
                       .AddOptionalNamedDoubleVector("FIBER2", 3)
                       .AddOptionalNamedDoubleVector("FIBER3", 3)
                       .AddOptionalNamedDouble("GROWTHTRIG")
                       .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoWeg6::SoWeg6(int id, int owner)
    : SoBase(id, owner), pstype_(INPAR::STR::PreStress::none), pstime_(0.0), time_(0.0)
{
  invJ_.resize(NUMGPT_WEG6);
  detJ_.resize(NUMGPT_WEG6);
  for (int i = 0; i < NUMGPT_WEG6; ++i)
  {
    detJ_[i] = 0.0;
    invJ_[i] = CORE::LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>(true);
  }

  Teuchos::RCP<const Teuchos::ParameterList> params =
      GLOBAL::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    pstype_ = PRESTRESS::GetType();
    pstime_ = PRESTRESS::GetPrestressTime();

    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        GLOBAL::Problem::Instance()->StructuralDynamicParams(), GetElementTypeString());
  }
  if (PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_WEG6, NUMGPT_WEG6));
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoWeg6::SoWeg6(const DRT::ELEMENTS::SoWeg6& old)
    : SoBase(old), detJ_(old.detJ_), pstype_(old.pstype_), pstime_(old.pstime_), time_(old.time_)
{
  invJ_.resize(old.invJ_.size());
  for (unsigned int i = 0; i < invJ_.size(); ++i)
  {
    invJ_[i] = old.invJ_[i];
  }

  if (PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoWeg6::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::SoWeg6(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::SoWeg6::Shape() const { return CORE::FE::CellType::wedge6; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoWeg6::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  SoBase::Pack(data);

  // Pack prestress
  AddtoPack(data, static_cast<int>(pstype_));
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (PRESTRESS::IsMulf(pstype_))
  {
    CORE::COMM::ParObject::AddtoPack(data, *prestress_);
  }

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const unsigned int size = invJ_.size();
  AddtoPack(data, size);
  for (unsigned int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoWeg6::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  SoBase::Unpack(basedata);
  // prestress_
  pstype_ = static_cast<INPAR::STR::PreStress>(ExtractInt(position, data));
  ExtractfromPack(position, data, pstime_);
  ExtractfromPack(position, data, time_);
  if (PRESTRESS::IsMulf(pstype_))
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_WEG6, NUMGPT_WEG6));
    prestress_->Unpack(tmpprestress);
  }

  // detJ_
  ExtractfromPack(position, data, detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size);
  for (int i = 0; i < size; ++i)
  {
    invJ_[i] = CORE::LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>(true);
    ExtractfromPack(position, data, invJ_[i]);
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoWeg6::Print(std::ostream& os) const
{
  os << "So_weg6 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}


std::vector<double> DRT::ELEMENTS::SoWeg6::ElementCenterRefeCoords()
{
  // update element geometry
  DRT::Node** nodes = Nodes();
  CORE::LINALG::Matrix<NUMNOD_WEG6, NUMDIM_WEG6> xrefe;  // material coord. of element
  for (int i = 0; i < NUMNOD_WEG6; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  const CORE::FE::CellType distype = Shape();
  CORE::LINALG::Matrix<NUMNOD_WEG6, 1> funct;
  // Element midpoint at r=s=1/3, t=0.0
  CORE::FE::shape_function_3D(funct, 1.0 / 3.0, 1.0 / 3.0, 0.0, distype);
  CORE::LINALG::Matrix<1, NUMDIM_WEG6> midpoint;
  // midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0, 0);
  centercoords[1] = midpoint(0, 1);
  centercoords[2] = midpoint(0, 2);
  return centercoords;
}
/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                maf 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoWeg6::VisNames(std::map<std::string, int>& names)
{
  SolidMaterial()->VisNames(names);

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         maf 07/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::SoWeg6::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, NUMGPT_WEG6, this->Id());
}


/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                             maf 04/07|
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::SoWeg6::Surfaces()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralSurface, DRT::Element>(
      CORE::COMM::buildSurfaces, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::SoWeg6::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralLine, DRT::Element>(
      CORE::COMM::buildLines, *this);
}

FOUR_C_NAMESPACE_CLOSE
