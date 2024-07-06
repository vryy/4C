/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element Type with ScaTra coupling

*----------------------------------------------------------------------*/

#include "4C_membrane_scatra_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_membrane_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneScatraTri3Type Discret::ELEMENTS::MembraneScatraTri3Type::instance_;

Discret::ELEMENTS::MembraneScatraTri3Type& Discret::ELEMENTS::MembraneScatraTri3Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::MembraneScatraTri3Type::create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri3>* object =
      new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri3>(-1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneScatraTri3Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA3" && eledistype == "TRI3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneScatraTri3Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri3>(id, owner));
  return ele;
}

void Discret::ELEMENTS::MembraneScatraTri3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_membrane;
  MembraneTri3Type::setup_element_definition(definitions_membrane);

  std::map<std::string, Input::LineDefinition>& defs_membrane = definitions_membrane["MEMBRANE3"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["MEMBRANESCATRA3"];

  defs["TRI3"] =
      Input::LineDefinition::Builder(defs_membrane["TRI3"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  TRI 6 Element                                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneScatraTri6Type Discret::ELEMENTS::MembraneScatraTri6Type::instance_;

Discret::ELEMENTS::MembraneScatraTri6Type& Discret::ELEMENTS::MembraneScatraTri6Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::MembraneScatraTri6Type::create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri6>* object =
      new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri6>(-1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneScatraTri6Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA6" && eledistype == "TRI6")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri6>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneScatraTri6Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri6>(id, owner));
  return ele;
}

void Discret::ELEMENTS::MembraneScatraTri6Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_membrane;
  MembraneTri6Type::setup_element_definition(definitions_membrane);

  std::map<std::string, Input::LineDefinition>& defs_membrane = definitions_membrane["MEMBRANE6"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["MEMBRANESCATRA6"];

  defs["TRI6"] =
      Input::LineDefinition::Builder(defs_membrane["TRI6"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneScatraQuad4Type Discret::ELEMENTS::MembraneScatraQuad4Type::instance_;

Discret::ELEMENTS::MembraneScatraQuad4Type& Discret::ELEMENTS::MembraneScatraQuad4Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::MembraneScatraQuad4Type::create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad4>* object =
      new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad4>(-1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneScatraQuad4Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA4" && eledistype == "QUAD4")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneScatraQuad4Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad4>(id, owner));
  return ele;
}

void Discret::ELEMENTS::MembraneScatraQuad4Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_membrane;
  MembraneQuad4Type::setup_element_definition(definitions_membrane);

  std::map<std::string, Input::LineDefinition>& defs_membrane = definitions_membrane["MEMBRANE4"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["MEMBRANESCATRA4"];

  defs["QUAD4"] =
      Input::LineDefinition::Builder(defs_membrane["QUAD4"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneScatraQuad9Type Discret::ELEMENTS::MembraneScatraQuad9Type::instance_;

Discret::ELEMENTS::MembraneScatraQuad9Type& Discret::ELEMENTS::MembraneScatraQuad9Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::MembraneScatraQuad9Type::create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad9>* object =
      new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad9>(-1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneScatraQuad9Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA9" && eledistype == "QUAD9")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneScatraQuad9Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad9>(id, owner));
  return ele;
}

void Discret::ELEMENTS::MembraneScatraQuad9Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_membrane;
  MembraneQuad9Type::setup_element_definition(definitions_membrane);

  std::map<std::string, Input::LineDefinition>& defs_membrane = definitions_membrane["MEMBRANE9"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["MEMBRANESCATRA9"];

  defs["QUAD9"] =
      Input::LineDefinition::Builder(defs_membrane["QUAD9"]).add_named_string("TYPE").build();
}

FOUR_C_NAMESPACE_CLOSE
