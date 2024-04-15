/*----------------------------------------------------------------------*/
/*! \file

\brief implementation of the quadratic NURBS 27 element

\level 2


*----------------------------------------------------------------------*/

#include "baci_so3_nurbs27.hpp"

#include "baci_comm_utils_factory.hpp"
#include "baci_global_data.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_lib_discret.hpp"
#include "baci_so3_line.hpp"
#include "baci_so3_nullspace.hpp"
#include "baci_so3_surface.hpp"
#include "baci_so3_utils.hpp"
#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::NURBS::SoNurbs27Type DRT::ELEMENTS::NURBS::SoNurbs27Type::instance_;

DRT::ELEMENTS::NURBS::SoNurbs27Type& DRT::ELEMENTS::NURBS::SoNurbs27Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::NURBS::SoNurbs27Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::NURBS::SoNurbs27(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::SoNurbs27Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::NURBS::SoNurbs27(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NURBS::SoNurbs27Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::NURBS::SoNurbs27(id, owner));
  return ele;
}


void DRT::ELEMENTS::NURBS::SoNurbs27Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::NURBS::SoNurbs27Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::NURBS::SoNurbs27Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["NURBS27"] = INPUT::LineDefinition::Builder()
                        .AddIntVector("NURBS27", 27)
                        .AddNamedInt("MAT")
                        .AddNamedIntVector("GP", 3)
                        .Build();
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::SoNurbs27::SoNurbs27(int id, int owner) : SoBase(id, owner)
{
  invJ_.resize(NUMGPT_SONURBS27, CORE::LINALG::Matrix<NUMDIM_SONURBS27, NUMDIM_SONURBS27>(true));
  detJ_.resize(NUMGPT_SONURBS27, 0.0);
  SetNurbsElement() = true;

  Teuchos::RCP<const Teuchos::ParameterList> params =
      GLOBAL::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        GLOBAL::Problem::Instance()->StructuralDynamicParams(), GetElementTypeString());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NURBS::SoNurbs27::SoNurbs27(const DRT::ELEMENTS::NURBS::SoNurbs27& old)
    : SoBase(old), detJ_(old.detJ_)
{
  invJ_.resize(old.invJ_.size());
  for (int i = 0; i < (int)invJ_.size(); ++i)
  {
    invJ_[i] = old.invJ_[i];
  }
  SetNurbsElement() = true;

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::NURBS::SoNurbs27::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::NURBS::SoNurbs27(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::NURBS::SoNurbs27::Shape() const
{
  return CORE::FE::CellType::nurbs27;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::SoNurbs27::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  SoBase::Pack(data);

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::SoNurbs27::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  SoBase::Unpack(basedata);
  // detJ_
  ExtractfromPack(position, data, detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, CORE::LINALG::Matrix<NUMDIM_SONURBS27, NUMDIM_SONURBS27>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::SoNurbs27::Print(std::ostream& os) const
{
  os << "So_nurbs27 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                                      |
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NURBS::SoNurbs27::Surfaces()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralSurface, SoNurbs27>(
      CORE::COMM::buildSurfaces, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                                        |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NURBS::SoNurbs27::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralLine, SoNurbs27>(
      CORE::COMM::buildLines, *this);
}

FOUR_C_NAMESPACE_CLOSE
