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

using namespace Core::FE;

Discret::ELEMENTS::RedInterAcinarDepType Discret::ELEMENTS::RedInterAcinarDepType::instance_;

Discret::ELEMENTS::RedInterAcinarDepType& Discret::ELEMENTS::RedInterAcinarDepType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::RedInterAcinarDepType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::RedInterAcinarDep* object = new Discret::ELEMENTS::RedInterAcinarDep(-1, -1);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::RedInterAcinarDepType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RED_ACINAR_INTER_DEP")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::RedInterAcinarDep(id, owner));
    return ele;
  }
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::RedInterAcinarDepType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::RedInterAcinarDep(id, owner));
  return ele;
}


/*----------------------------------------------------------------------*
 |  setup_element_definition                                              |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedInterAcinarDepType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["RED_ACINAR_INTER_DEP"];

  defs["LINE2"] =
      Input::LineDefinition::Builder().AddIntVector("LINE2", 2).AddNamedInt("MAT").Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::RedInterAcinarDep::RedInterAcinarDep(int id, int owner)
    : Core::Elements::Element(id, owner)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::RedInterAcinarDep::RedInterAcinarDep(
    const Discret::ELEMENTS::RedInterAcinarDep& old)
    : Core::Elements::Element(old), elem_params_(old.elem_params_), generation_(old.generation_)
{
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedInterAcinarDep and return pointer     |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::RedInterAcinarDep::Clone() const
{
  Discret::ELEMENTS::RedInterAcinarDep* newelement =
      new Discret::ELEMENTS::RedInterAcinarDep(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::RedInterAcinarDep::Shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      return Core::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedInterAcinarDep::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // add base class Element
  Element::Pack(data);

  std::map<std::string, double>::const_iterator it;

  add_to_pack(data, (int)(elem_params_.size()));
  for (it = elem_params_.begin(); it != elem_params_.end(); it++)
  {
    add_to_pack(data, it->first);
    add_to_pack(data, it->second);
  }

  add_to_pack(data, generation_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedInterAcinarDep::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::Unpack(basedata);

  std::map<std::string, double> it;
  int n = 0;

  extract_from_pack(position, data, n);

  for (int i = 0; i < n; i++)
  {
    std::string name;
    double val;
    extract_from_pack(position, data, name);
    extract_from_pack(position, data, val);
    elem_params_[name] = val;
  }

  // extract generation
  extract_from_pack(position, data, generation_);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}



/*----------------------------------------------------------------------*
 |  Print this element (public)                             ismail 01/10|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedInterAcinarDep::Print(std::ostream& os) const
{
  os << "RedInterAcinarDep ";
  Element::Print(os);

  return;
}


/*----------------------------------------------------------------------*
 |  Return names of visualization data                     ismail 01/10 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedInterAcinarDep::VisNames(std::map<std::string, int>& names) { return; }

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::RedInterAcinarDep::VisData(
    const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::VisData(name, data)) return true;

  return false;
}


/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 04/10 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedInterAcinarDep::getParams(std::string name, double& var)
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
void Discret::ELEMENTS::RedInterAcinarDep::getParams(std::string name, int& var)
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
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::RedInterAcinarDep::Lines()
{
  FOUR_C_ASSERT(NumLine() == 1, "RED_AIRWAY element must have one and only one line");

  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
