/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_bele_vele3.hpp"

#include "baci_comm_utils_factory.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_lib_discret.hpp"
#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

DRT::ELEMENTS::Vele3Type DRT::ELEMENTS::Vele3Type::instance_;

DRT::ELEMENTS::Vele3Type& DRT::ELEMENTS::Vele3Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::Vele3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Vele3* object = new DRT::ELEMENTS::Vele3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Vele3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "VELE3")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Vele3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Vele3Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Vele3(id, owner));
  return ele;
}


void DRT::ELEMENTS::Vele3Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::Vele3Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  dserror("method ComputeNullSpace not implemented for element type vele3!");
  return nullspace;
}

void DRT::ELEMENTS::Vele3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["VELE3"];

  defs["HEX8"] = INPUT::LineDefinition::Builder().AddIntVector("HEX8", 8).Build();
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Vele3SurfaceType::Create(const int id, const int owner)
{
  // return Teuchos::rcp( new Vele3Surface( id, owner ) );
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Vele3LineType::Create(const int id, const int owner)
{
  // return Teuchos::rcp( new Vele3Line( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3::Vele3(int id, int owner) : DRT::Element(id, owner) { return; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3::Vele3(const DRT::ELEMENTS::Vele3& old) : DRT::Element(old) { return; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Vele3::Clone() const
{
  DRT::ELEMENTS::Vele3* newelement = new DRT::ELEMENTS::Vele3(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Vele3::Shape() const
{
  switch (NumNode())
  {
    case 4:
      return CORE::FE::CellType::tet4;
    case 5:
      return CORE::FE::CellType::pyramid5;
    case 6:
      return CORE::FE::CellType::wedge6;
    case 8:
      return CORE::FE::CellType::hex8;
    case 10:
      return CORE::FE::CellType::tet10;
    case 15:
      return CORE::FE::CellType::wedge15;
    case 20:
      return CORE::FE::CellType::hex20;
    case 27:
      return CORE::FE::CellType::hex27;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3::Print(std::ostream& os) const
{
  os << "Vele3 " << CORE::FE::CellTypeToString(Shape());
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Vele3::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<Vele3Line, Vele3>(CORE::COMM::buildLines, *this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Vele3::Surfaces()
{
  return CORE::COMM::ElementBoundaryFactory<Vele3Surface, Vele3>(CORE::COMM::buildSurfaces, *this);
}



/*----------------------------------------------------------------------*
 |  get optimal gauss rule (public)                          u.may 05/09|
 *----------------------------------------------------------------------*/
CORE::FE::GaussRule3D DRT::ELEMENTS::Vele3::getOptimalGaussrule(
    const CORE::FE::CellType& distype) const
{
  CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::undefined;
  switch (distype)
  {
    case CORE::FE::CellType::hex8:
      rule = CORE::FE::GaussRule3D::hex_8point;
      break;
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
      rule = CORE::FE::GaussRule3D::hex_27point;
      break;
    case CORE::FE::CellType::tet4:
      rule = CORE::FE::GaussRule3D::tet_4point;
      break;
    case CORE::FE::CellType::tet10:
      rule = CORE::FE::GaussRule3D::tet_10point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Vele3::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  return true;
}

FOUR_C_NAMESPACE_CLOSE
