/*----------------------------------------------------------------------*/
/*!

\brief 2D boundary elment

\maintainer Amadeus Gebauer

\level 2
*/
/*----------------------------------------------------------------------*/

#include "bele2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_utils_nullspace.H"


DRT::ELEMENTS::Bele2Type DRT::ELEMENTS::Bele2Type::instance_;

DRT::ELEMENTS::Bele2Type& DRT::ELEMENTS::Bele2Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Bele2Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Bele2* object = new DRT::ELEMENTS::Bele2(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Bele2Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "BELE2")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Bele2(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Bele2Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Bele2(id, owner));
  return ele;
}


void DRT::ELEMENTS::Bele2Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 2;
  dimns = 2;
  nv = 2;
}

void DRT::ELEMENTS::Bele2Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure2DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Bele2Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["BELE2"];

  defs["LINE2"].AddIntVector("LINE2", 2);

  defs["LINE3"].AddIntVector("LINE3", 3);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2::Bele2(int id, int owner) : DRT::Element(id, owner) { return; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2::Bele2(const DRT::ELEMENTS::Bele2& old) : DRT::Element(old) { return; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Bele2::Clone() const
{
  DRT::ELEMENTS::Bele2* newelement = new DRT::ELEMENTS::Bele2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Bele2::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return line2;
    case 3:
      return line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele2::Pack(DRT::PackBuffer& data) const
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
void DRT::ELEMENTS::Bele2::Unpack(const std::vector<char>& data)
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
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2::~Bele2() { return; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele2::Print(std::ostream& os) const
{
  os << "Bele2 " << DRT::DistypeToString(Shape());
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Bele2::Lines()
{
  std::vector<Teuchos::RCP<DRT::Element>> lines(1);
  lines[0] = Teuchos::rcp(this, false);
  return lines;
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Bele2::Surfaces()
{
  std::vector<Teuchos::RCP<DRT::Element>> surfaces(0);
  return surfaces;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Bele2::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  return true;
}
