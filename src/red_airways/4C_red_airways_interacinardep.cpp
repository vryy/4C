/*---------------------------------------------------------------------*/
/*! \file

\brief Implements an inter-acinar dependency (element)


\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

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
Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::RedInterAcinarDepType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RED_ACINAR_INTER_DEP")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::RedInterAcinarDep(id, owner));
    return ele;
  }
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::RedInterAcinarDepType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::RedInterAcinarDep(id, owner));
  return ele;
}


/*----------------------------------------------------------------------*
 |  setup_element_definition                                              |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedInterAcinarDepType::setup_element_definition(
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
DRT::ELEMENTS::RedInterAcinarDep::RedInterAcinarDep(int id, int owner)
    : CORE::Elements::Element(id, owner)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedInterAcinarDep::RedInterAcinarDep(const DRT::ELEMENTS::RedInterAcinarDep& old)
    : CORE::Elements::Element(old), elem_params_(old.elem_params_), generation_(old.generation_)
{
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedInterAcinarDep and return pointer     |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::RedInterAcinarDep::Clone() const
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
  switch (num_node())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
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

  AddtoPack(data, (int)(elem_params_.size()));
  for (it = elem_params_.begin(); it != elem_params_.end(); it++)
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
    elem_params_[name] = val;
  }

  // extract generation
  ExtractfromPack(position, data, generation_);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

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
  if (CORE::Elements::Element::VisData(name, data)) return true;

  return false;
}


/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 04/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedInterAcinarDep::getParams(std::string name, double& var)
{
  std::map<std::string, double>::iterator it;
  it = elem_params_.find(name);
  if (it == elem_params_.end())
  {
    FOUR_C_THROW("[%s] is not found with in the element variables", name.c_str());
    exit(1);
  }
  var = elem_params_[name];
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
    FOUR_C_THROW("[%s] is not found with in the element INT variables", name.c_str());
    exit(1);
  }
}


/*----------------------------------------------------------------------*
 |  Get vector of lines (public)                           ismail  02/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::RedInterAcinarDep::Lines()
{
  FOUR_C_ASSERT(NumLine() == 1, "RED_AIRWAY element must have one and only one line");

  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
