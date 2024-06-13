/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element including scatra functionality.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1_poro_scatra_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_w1_poro_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::WallQuad4PoroScatraType Discret::ELEMENTS::WallQuad4PoroScatraType::instance_;

Discret::ELEMENTS::WallQuad4PoroScatraType& Discret::ELEMENTS::WallQuad4PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallQuad4PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad4>* object =
      new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad4PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ4POROSCATRA")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad4PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad4>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::WallQuad4PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  WallQuad4PoroType::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALLQ4PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLQ4POROSCATRA"];

  defs["QUAD4"] =
      Input::LineDefinition::Builder(defs_wall["QUAD4"]).add_named_string("TYPE").build();
}


/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::WallQuad9PoroScatraType Discret::ELEMENTS::WallQuad9PoroScatraType::instance_;

Discret::ELEMENTS::WallQuad9PoroScatraType& Discret::ELEMENTS::WallQuad9PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallQuad9PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad9>* object =
      new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad9PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ9POROSCATRA")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad9PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad9>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::WallQuad9PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  WallQuad9PoroType::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALLQ9PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLQ9POROSCATRA"];

  defs["QUAD9"] =
      Input::LineDefinition::Builder(defs_wall["QUAD9"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  NURBS 4 Element                                       schmidt 09/17 |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::WallNurbs4PoroScatraType Discret::ELEMENTS::WallNurbs4PoroScatraType::instance_;

Discret::ELEMENTS::WallNurbs4PoroScatraType& Discret::ELEMENTS::WallNurbs4PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallNurbs4PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs4>* object =
      new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs4>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallNurbs4PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLN4POROSCATRA")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallNurbs4PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs4>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::WallNurbs4PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  WallNurbs4PoroType::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALLN4PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLN4POROSCATRA"];

  defs["NURBS4"] =
      Input::LineDefinition::Builder(defs_wall["NURBS4"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  NURBS 9 Element                                       schmidt 09/17 |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::WallNurbs9PoroScatraType Discret::ELEMENTS::WallNurbs9PoroScatraType::instance_;

Discret::ELEMENTS::WallNurbs9PoroScatraType& Discret::ELEMENTS::WallNurbs9PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallNurbs9PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs9>* object =
      new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs9>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallNurbs9PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLN9POROSCATRA")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallNurbs9PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs9>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::WallNurbs9PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  WallNurbs9PoroType::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALLN9PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLN9POROSCATRA"];

  defs["NURBS9"] =
      Input::LineDefinition::Builder(defs_wall["NURBS9"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::WallTri3PoroScatraType Discret::ELEMENTS::WallTri3PoroScatraType::instance_;

Discret::ELEMENTS::WallTri3PoroScatraType& Discret::ELEMENTS::WallTri3PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallTri3PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::tri3>* object =
      new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallTri3PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLT3POROSCATRA")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallTri3PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::tri3>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::WallTri3PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  WallTri3PoroType::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALLT3PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLT3POROSCATRA"];

  defs["TRI3"] = Input::LineDefinition::Builder(defs_wall["TRI3"]).add_named_string("TYPE").build();
}

FOUR_C_NAMESPACE_CLOSE
