/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_bele_vele3.H"

#include "baci_io_linedefinition.H"
#include "baci_lib_discret.H"
#include "baci_lib_utils_factory.H"
#include "baci_utils_exceptions.H"

DRT::ELEMENTS::Vele3Type DRT::ELEMENTS::Vele3Type::instance_;

DRT::ELEMENTS::Vele3Type& DRT::ELEMENTS::Vele3Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Vele3Type::Create(const std::vector<char>& data)
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
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["VELE3"];

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
DRT::Element::DiscretizationType DRT::ELEMENTS::Vele3::Shape() const
{
  switch (NumNode())
  {
    case 4:
      return DiscretizationType::tet4;
    case 5:
      return DiscretizationType::pyramid5;
    case 6:
      return DiscretizationType::wedge6;
    case 8:
      return DiscretizationType::hex8;
    case 10:
      return DiscretizationType::tet10;
    case 15:
      return DiscretizationType::wedge15;
    case 20:
      return DiscretizationType::hex20;
    case 27:
      return DiscretizationType::hex27;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
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
  os << "Vele3 " << DRT::DistypeToString(Shape());
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Vele3::Lines()
{
  return DRT::UTILS::ElementBoundaryFactory<Vele3Line, Vele3>(DRT::UTILS::buildLines, *this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Vele3::Surfaces()
{
  return DRT::UTILS::ElementBoundaryFactory<Vele3Surface, Vele3>(DRT::UTILS::buildSurfaces, *this);
}



/*----------------------------------------------------------------------*
 |  get optimal gauss rule (public)                          u.may 05/09|
 *----------------------------------------------------------------------*/
CORE::DRT::UTILS::GaussRule3D DRT::ELEMENTS::Vele3::getOptimalGaussrule(
    const DRT::Element::DiscretizationType& distype) const
{
  CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
  switch (distype)
  {
    case DRT::Element::DiscretizationType::hex8:
      rule = CORE::DRT::UTILS::GaussRule3D::hex_8point;
      break;
    case DRT::Element::DiscretizationType::hex20:
    case DRT::Element::DiscretizationType::hex27:
      rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
      break;
    case DRT::Element::DiscretizationType::tet4:
      rule = CORE::DRT::UTILS::GaussRule3D::tet_4point;
      break;
    case DRT::Element::DiscretizationType::tet10:
      rule = CORE::DRT::UTILS::GaussRule3D::tet_10point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Vele3::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  return true;
}
