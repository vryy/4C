/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D wall element for solid-part of porous medium using p1 (mixed)
 approach including scatra functionality.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1_poro_p1_scatra_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_w1_poro_p1_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::WallQuad4PoroP1ScatraType
    Discret::ELEMENTS::WallQuad4PoroP1ScatraType::instance_;

Discret::ELEMENTS::WallQuad4PoroP1ScatraType&
Discret::ELEMENTS::WallQuad4PoroP1ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallQuad4PoroP1ScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad4>* object =
      new Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad4PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ4POROP1SCATRA")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad4PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad4>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::WallQuad4PoroP1ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wallporo;
  WallQuad4PoroP1Type::setup_element_definition(definitions_wallporo);

  std::map<std::string, Input::LineDefinition>& defs_wallporo =
      definitions_wallporo["WALLQ4POROP1"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLQ4POROP1SCATRA"];

  defs["QUAD4"] =
      Input::LineDefinition::Builder(defs_wallporo["QUAD4"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::WallQuad9PoroP1ScatraType
    Discret::ELEMENTS::WallQuad9PoroP1ScatraType::instance_;

Discret::ELEMENTS::WallQuad9PoroP1ScatraType&
Discret::ELEMENTS::WallQuad9PoroP1ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallQuad9PoroP1ScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad9>* object =
      new Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad9PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ9POROP1SCATRA")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad9PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad9>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::WallQuad9PoroP1ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wallporo;
  WallQuad9PoroP1Type::setup_element_definition(definitions_wallporo);

  std::map<std::string, Input::LineDefinition>& defs_wallporo =
      definitions_wallporo["WALLQ9POROP1"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLQ9POROP1SCATRA"];

  defs["QUAD9"] =
      Input::LineDefinition::Builder(defs_wallporo["QUAD9"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::WallTri3PoroP1ScatraType Discret::ELEMENTS::WallTri3PoroP1ScatraType::instance_;

Discret::ELEMENTS::WallTri3PoroP1ScatraType& Discret::ELEMENTS::WallTri3PoroP1ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallTri3PoroP1ScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::tri3>* object =
      new Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallTri3PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLT3POROP1SCATRA")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallTri3PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::tri3>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::WallTri3PoroP1ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wallporo;
  WallTri3PoroP1Type::setup_element_definition(definitions_wallporo);

  std::map<std::string, Input::LineDefinition>& defs_wallporo =
      definitions_wallporo["WALLT3POROP1"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLT3POROP1SCATRA"];

  defs["TRI3"] =
      Input::LineDefinition::Builder(defs_wallporo["TRI3"]).add_named_string("TYPE").build();
}

FOUR_C_NAMESPACE_CLOSE
