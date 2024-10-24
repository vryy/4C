// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_membrane_scatra_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_membrane_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
Discret::Elements::MembraneScatraTri3Type Discret::Elements::MembraneScatraTri3Type::instance_;

Discret::Elements::MembraneScatraTri3Type& Discret::Elements::MembraneScatraTri3Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::MembraneScatraTri3Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::MembraneScatra<Core::FE::CellType::tri3>* object =
      new Discret::Elements::MembraneScatra<Core::FE::CellType::tri3>(-1, -1);
  object->unpack(buffer);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneScatraTri3Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA3" && eledistype == "TRI3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::Elements::MembraneScatra<Core::FE::CellType::tri3>>(id, owner);
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneScatraTri3Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::Elements::MembraneScatra<Core::FE::CellType::tri3>>(id, owner);
  return ele;
}

void Discret::Elements::MembraneScatraTri3Type::setup_element_definition(
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
Discret::Elements::MembraneScatraTri6Type Discret::Elements::MembraneScatraTri6Type::instance_;

Discret::Elements::MembraneScatraTri6Type& Discret::Elements::MembraneScatraTri6Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::MembraneScatraTri6Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::MembraneScatra<Core::FE::CellType::tri6>* object =
      new Discret::Elements::MembraneScatra<Core::FE::CellType::tri6>(-1, -1);
  object->unpack(buffer);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneScatraTri6Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA6" && eledistype == "TRI6")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::Elements::MembraneScatra<Core::FE::CellType::tri6>>(id, owner);
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneScatraTri6Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::Elements::MembraneScatra<Core::FE::CellType::tri6>>(id, owner);
  return ele;
}

void Discret::Elements::MembraneScatraTri6Type::setup_element_definition(
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
Discret::Elements::MembraneScatraQuad4Type Discret::Elements::MembraneScatraQuad4Type::instance_;

Discret::Elements::MembraneScatraQuad4Type& Discret::Elements::MembraneScatraQuad4Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::MembraneScatraQuad4Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::MembraneScatra<Core::FE::CellType::quad4>* object =
      new Discret::Elements::MembraneScatra<Core::FE::CellType::quad4>(-1, -1);
  object->unpack(buffer);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneScatraQuad4Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA4" && eledistype == "QUAD4")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::Elements::MembraneScatra<Core::FE::CellType::quad4>>(id, owner);
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneScatraQuad4Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::Elements::MembraneScatra<Core::FE::CellType::quad4>>(id, owner);
  return ele;
}

void Discret::Elements::MembraneScatraQuad4Type::setup_element_definition(
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
Discret::Elements::MembraneScatraQuad9Type Discret::Elements::MembraneScatraQuad9Type::instance_;

Discret::Elements::MembraneScatraQuad9Type& Discret::Elements::MembraneScatraQuad9Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::MembraneScatraQuad9Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::MembraneScatra<Core::FE::CellType::quad9>* object =
      new Discret::Elements::MembraneScatra<Core::FE::CellType::quad9>(-1, -1);
  object->unpack(buffer);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneScatraQuad9Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA9" && eledistype == "QUAD9")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::Elements::MembraneScatra<Core::FE::CellType::quad9>>(id, owner);
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneScatraQuad9Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::Elements::MembraneScatra<Core::FE::CellType::quad9>>(id, owner);
  return ele;
}

void Discret::Elements::MembraneScatraQuad9Type::setup_element_definition(
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
