/*---------------------------------------------------------------------*/
/*! \file

\brief Incomplete! - Purpose: Implements RedAirBloodScatraLine3 element


\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

using namespace Core::FE;

Discret::ELEMENTS::RedAirBloodScatraLine3Type
    Discret::ELEMENTS::RedAirBloodScatraLine3Type::instance_;

Discret::ELEMENTS::RedAirBloodScatraLine3Type&
Discret::ELEMENTS::RedAirBloodScatraLine3Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::RedAirBloodScatraLine3Type::create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::RedAirBloodScatraLine3* object =
      new Discret::ELEMENTS::RedAirBloodScatraLine3(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::RedAirBloodScatraLine3Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RED_AIR_BLOOD_SCATRA_LINE3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::RedAirBloodScatraLine3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::RedAirBloodScatraLine3Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::RedAirBloodScatraLine3(id, owner));
  return ele;
}


void Discret::ELEMENTS::RedAirBloodScatraLine3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["RED_AIR_BLOOD_SCATRA_LINE3"];

  defs["LINE3"] = Input::LineDefinition::Builder()
                      .add_int_vector("LINE3", 3)
                      .add_named_double("DiffusionCoefficient")
                      .add_named_double("WallThickness")
                      .add_named_double("PercentageOfDiffusionArea")
                      .build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 05/13|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::RedAirBloodScatraLine3::RedAirBloodScatraLine3(int id, int owner)
    : Core::Elements::Element(id, owner)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 05/13|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::RedAirBloodScatraLine3::RedAirBloodScatraLine3(
    const Discret::ELEMENTS::RedAirBloodScatraLine3& old)
    : Core::Elements::Element(old), elem_params_(old.elem_params_), generation_(old.generation_)
{
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedAirBloodScatraLine3 and return pointer             |
 |  to it                                                      (public) |
 |                                                         ismail 05/13 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::RedAirBloodScatraLine3::clone() const
{
  Discret::ELEMENTS::RedAirBloodScatraLine3* newelement =
      new Discret::ELEMENTS::RedAirBloodScatraLine3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 05/13 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::RedAirBloodScatraLine3::shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      return Core::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 05/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedAirBloodScatraLine3::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Element
  Element::pack(data);


  std::map<std::string, double>::const_iterator it;

  add_to_pack(data, (int)(elem_params_.size()));
  for (it = elem_params_.begin(); it != elem_params_.end(); it++)
  {
    add_to_pack(data, it->first);
    add_to_pack(data, it->second);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 05/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedAirBloodScatraLine3::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);

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


  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 05/13|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedAirBloodScatraLine3::print(std::ostream& os) const
{
  os << "RedAirBloodScatraLine3 ";
  Element::print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data                     ismail 05/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedAirBloodScatraLine3::vis_names(std::map<std::string, int>& names)
{
  // Put the owner of this element into the file (use base class method for this)
  Core::Elements::Element::vis_names(names);
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::RedAirBloodScatraLine3::vis_data(
    const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return false;
}



/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 04/10 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RedAirBloodScatraLine3::get_params(std::string name, double& var)
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
void Discret::ELEMENTS::RedAirBloodScatraLine3::get_params(std::string name, int& var)
{
  //  if (name == "Generation")
  //  {
  //    var = generation_;
  //  }
  //  else
  {
    FOUR_C_THROW("[%s] is not found with in the element INT variables", name.c_str());
    exit(1);
  }
}

/*----------------------------------------------------------------------*
 |  get vector of lines              (public)              ismail  02/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::RedAirBloodScatraLine3::lines()
{
  FOUR_C_ASSERT(num_line() == 1, "RED_AIRWAY element must have one and only one line");

  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
