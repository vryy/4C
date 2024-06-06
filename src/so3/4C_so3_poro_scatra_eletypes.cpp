/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the 3D solid-poro element including scatra functionality

 \level 2

*----------------------------------------------------------------------*/

#include "4C_so3_poro_scatra_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_poro_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8PoroScatraType Discret::ELEMENTS::SoHex8PoroScatraType::instance_;

Discret::ELEMENTS::SoHex8PoroScatraType& Discret::ELEMENTS::SoHex8PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoHex8PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>* object =
      new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoHex8PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex8;
  SoHex8PoroType::setup_element_definition(definitions_hex8);

  std::map<std::string, Input::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = Input::LineDefinition::Builder(defs_hex8["HEX8"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/


Discret::ELEMENTS::SoTet4PoroScatraType Discret::ELEMENTS::SoTet4PoroScatraType::instance_;

Discret::ELEMENTS::SoTet4PoroScatraType& Discret::ELEMENTS::SoTet4PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoTet4PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>* object =
      new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoTet4PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_tet4;
  SoTet4PoroType::setup_element_definition(definitions_tet4);

  std::map<std::string, Input::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = Input::LineDefinition::Builder(defs_tet4["TET4"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  HEX 27 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/


Discret::ELEMENTS::SoHex27PoroScatraType Discret::ELEMENTS::SoHex27PoroScatraType::instance_;

Discret::ELEMENTS::SoHex27PoroScatraType& Discret::ELEMENTS::SoHex27PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoHex27PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>* object =
      new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex27PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex27PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoHex27PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex27;
  SoHex27PoroType::setup_element_definition(definitions_hex27);

  std::map<std::string, Input::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX27"] =
      Input::LineDefinition::Builder(defs_hex27["HEX27"]).AddNamedString("TYPE").Build();
}


/*----------------------------------------------------------------------*
 |  TET 10 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/


Discret::ELEMENTS::SoTet10PoroScatraType Discret::ELEMENTS::SoTet10PoroScatraType::instance_;

Discret::ELEMENTS::SoTet10PoroScatraType& Discret::ELEMENTS::SoTet10PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoTet10PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>* object =
      new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet10PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet10PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoTet10PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_tet10;
  SoTet10PoroType::setup_element_definition(definitions_tet10);

  std::map<std::string, Input::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET10"] =
      Input::LineDefinition::Builder(defs_tet10["TET10"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  NURBS 27 Element                                      schmidt 09/17 |
 *----------------------------------------------------------------------*/


Discret::ELEMENTS::SoNurbs27PoroScatraType Discret::ELEMENTS::SoNurbs27PoroScatraType::instance_;

Discret::ELEMENTS::SoNurbs27PoroScatraType& Discret::ELEMENTS::SoNurbs27PoroScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoNurbs27PoroScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::Nurbs::SoNurbs27,
      Core::FE::CellType::nurbs27>* object =
      new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::Nurbs::SoNurbs27,
          Core::FE::CellType::nurbs27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoNurbs27PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::Nurbs::SoNurbs27,
            Core::FE::CellType::nurbs27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoNurbs27PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::So3PoroScatra<Discret::ELEMENTS::Nurbs::SoNurbs27,
          Core::FE::CellType::nurbs27>(id, owner));
  return ele;
}

void Discret::ELEMENTS::SoNurbs27PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_nurbs27;
  SoNurbs27PoroType::setup_element_definition(definitions_nurbs27);

  std::map<std::string, Input::LineDefinition>& defs_nurbs27 = definitions_nurbs27["SONURBS27PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["NURBS27"] =
      Input::LineDefinition::Builder(defs_nurbs27["NURBS27"]).AddNamedString("TYPE").Build();
}

FOUR_C_NAMESPACE_CLOSE
