/*----------------------------------------------------------------------------*/
/*! \file
\brief a 2D solid-wall element with ScaTra coupling.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1_scatra.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Wall1ScatraType Discret::ELEMENTS::Wall1ScatraType::instance_;

Discret::ELEMENTS::Wall1ScatraType& Discret::ELEMENTS::Wall1ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::Wall1ScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1Scatra* object = new Discret::ELEMENTS::Wall1Scatra(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Wall1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLSCATRA")
  {
    if (eledistype != "NURBS4" and eledistype != "NURBS9")
    {
      return Teuchos::rcp(new Discret::ELEMENTS::Wall1Scatra(id, owner));
    }
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Wall1ScatraType::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::Wall1Scatra(id, owner));
}

void Discret::ELEMENTS::Wall1ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLSCATRA"];

  for (const auto& [key, wall_line_def] : defs_wall)
  {
    defs[key] = Input::LineDefinition::Builder(wall_line_def).add_named_string("TYPE").build();
  }
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 01/14/|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Wall1Scatra::Wall1Scatra(int id, int owner)
    : Wall1(id, owner), impltype_(Inpar::ScaTra::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 01/14|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Wall1Scatra::Wall1Scatra(const Discret::ELEMENTS::Wall1Scatra& old)
    : Wall1(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                            vuong 01/14 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Wall1Scatra::Clone() const
{
  Discret::ELEMENTS::Wall1Scatra* newelement = new Discret::ELEMENTS::Wall1Scatra(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            vuong 01/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Wall1Scatra::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // pack scalar transport impltype
  add_to_pack(data, impltype_);

  // add base class Element
  Wall1::pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            vuong 01/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Wall1Scatra::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract scalar transport impltype
  impltype_ = static_cast<Inpar::ScaTra::ImplType>(extract_int(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Wall1::unpack(basedata);
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 01/14|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Wall1Scatra::Print(std::ostream& os) const
{
  os << "Wall1_Scatra ";
  Wall1::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Wall1Scatra::ReadElement(
    const std::string& eletype, const std::string& eledistype, Input::LineDefinition* linedef)
{
  // read base element
  Wall1::ReadElement(eletype, eledistype, linedef);

  // read scalar transport implementation type
  std::string impltype;
  linedef->extract_string("TYPE", impltype);

  if (impltype == "Undefined")
    impltype_ = Inpar::ScaTra::impltype_undefined;
  else if (impltype == "AdvReac")
    impltype_ = Inpar::ScaTra::impltype_advreac;
  else if (impltype == "CardMono")
    impltype_ = Inpar::ScaTra::impltype_cardiac_monodomain;
  else if (impltype == "Chemo")
    impltype_ = Inpar::ScaTra::impltype_chemo;
  else if (impltype == "ChemoReac")
    impltype_ = Inpar::ScaTra::impltype_chemoreac;
  else if (impltype == "Loma")
    impltype_ = Inpar::ScaTra::impltype_loma;
  else if (impltype == "RefConcReac")
    impltype_ = Inpar::ScaTra::impltype_refconcreac;
  else if (impltype == "Std")
    impltype_ = Inpar::ScaTra::impltype_std;
  else
    FOUR_C_THROW("Invalid implementation type for Wall1_Scatra elements!");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
