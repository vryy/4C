/*---------------------------------------------------------------------*/
/*! \file

\brief Implements an inter-acinar dependency (element)


\level 3

*/
/*---------------------------------------------------------------------*/

#include "baci_io_linedefinition.hpp"
#include "baci_lib_discret.hpp"
#include "baci_red_airways_elementbase.hpp"
#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN

using namespace CORE::FE;

DRT::ELEMENTS::RedInterAcinarDepType DRT::ELEMENTS::RedInterAcinarDepType::instance_;

DRT::ELEMENTS::RedInterAcinarDepType& DRT::ELEMENTS::RedInterAcinarDepType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::RedInterAcinarDepType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::RedInterAcinarDep* object = new DRT::ELEMENTS::RedInterAcinarDep(-1, -1);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedInterAcinarDepType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RED_ACINAR_INTER_DEP")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::RedInterAcinarDep(id, owner));
    return ele;
  }
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedInterAcinarDepType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::RedInterAcinarDep(id, owner));
  return ele;
}


/*----------------------------------------------------------------------*
 |  SetupElementDefinition                                              |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedInterAcinarDepType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["RED_ACINAR_INTER_DEP"];

  defs["LINE2"] =
      INPUT::LineDefinition::Builder().AddIntVector("LINE2", 2).AddNamedInt("MAT").Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedInterAcinarDep::RedInterAcinarDep(int id, int owner) : DRT::Element(id, owner) {}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedInterAcinarDep::RedInterAcinarDep(const DRT::ELEMENTS::RedInterAcinarDep& old)
    : DRT::Element(old), elemParams_(old.elemParams_), generation_(old.generation_)
{
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedInterAcinarDep and return pointer     |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::RedInterAcinarDep::Clone() const
{
  DRT::ELEMENTS::RedInterAcinarDep* newelement = new DRT::ELEMENTS::RedInterAcinarDep(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::RedInterAcinarDep::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
      break;
  }
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedInterAcinarDep::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  Element::Pack(data);

  std::map<std::string, double>::const_iterator it;

  AddtoPack(data, (int)(elemParams_.size()));
  for (it = elemParams_.begin(); it != elemParams_.end(); it++)
  {
    AddtoPack(data, it->first);
    AddtoPack(data, it->second);
  }

  AddtoPack(data, generation_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedInterAcinarDep::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  std::map<std::string, double> it;
  int n = 0;

  ExtractfromPack(position, data, n);

  for (int i = 0; i < n; i++)
  {
    std::string name;
    double val;
    ExtractfromPack(position, data, name);
    ExtractfromPack(position, data, val);
    elemParams_[name] = val;
  }

  // extract generation
  ExtractfromPack(position, data, generation_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}



/*----------------------------------------------------------------------*
 |  Print this element (public)                             ismail 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedInterAcinarDep::Print(std::ostream& os) const
{
  os << "RedInterAcinarDep ";
  Element::Print(os);

  return;
}


/*----------------------------------------------------------------------*
 |  Return names of visualization data                     ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedInterAcinarDep::VisNames(std::map<std::string, int>& names) { return; }

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedInterAcinarDep::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return false;
}


/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 04/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedInterAcinarDep::getParams(std::string name, double& var)
{
  std::map<std::string, double>::iterator it;
  it = elemParams_.find(name);
  if (it == elemParams_.end())
  {
    dserror("[%s] is not found with in the element variables", name.c_str());
    exit(1);
  }
  var = elemParams_[name];
}


/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 03/11 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedInterAcinarDep::getParams(std::string name, int& var)
{
  if (name == "Generation")
  {
    var = generation_;
  }
  else
  {
    dserror("[%s] is not found with in the element INT variables", name.c_str());
    exit(1);
  }
}


/*----------------------------------------------------------------------*
 |  Get vector of lines (public)                           ismail  02/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::RedInterAcinarDep::Lines()
{
  dsassert(NumLine() == 1, "RED_AIRWAY element must have one and only one line");

  return {Teuchos::rcpFromRef(*this)};
}

BACI_NAMESPACE_CLOSE
