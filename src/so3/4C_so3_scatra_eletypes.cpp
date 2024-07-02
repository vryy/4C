/*----------------------------------------------------------------------*/
/*! \file

\brief Solid-scatra elements types

\level 2


*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_so3_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8ScatraType Discret::ELEMENTS::SoHex8ScatraType::instance_;

Discret::ELEMENTS::SoHex8ScatraType& Discret::ELEMENTS::SoHex8ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoHex8ScatraType::Create(
    const std::vector<char>& data)
{
  auto* object =
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(-1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoHex8ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex8;
  SoHex8Type::setup_element_definition(definitions_hex8);

  std::map<std::string, Input::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8_DEPRECATED"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = Input::LineDefinition::Builder(defs_hex8["HEX8"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex8ScatraType::initialize(Core::FE::Discretization& dis)
{
  SoHex8Type::initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8_scatra* failed");
    actele->InitElement();
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  HEX 8 fbar Element                                       Thon 12/14 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8fbarScatraType Discret::ELEMENTS::SoHex8fbarScatraType::instance_;

Discret::ELEMENTS::SoHex8fbarScatraType& Discret::ELEMENTS::SoHex8fbarScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoHex8fbarScatraType::Create(
    const std::vector<char>& data)
{
  auto* object =
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>(
          -1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8fbarScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8fbarScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoHex8fbarScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex8;
  SoHex8fbarType::setup_element_definition(definitions_hex8);

  std::map<std::string, Input::LineDefinition>& defs_hex8 =
      definitions_hex8["SOLIDH8FBAR_DEPRECATED"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = Input::LineDefinition::Builder(defs_hex8["HEX8"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex8fbarScatraType::initialize(Core::FE::Discretization& dis)
{
  SoHex8fbarType::initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8fbar_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  HEX 27 Solid Scatra Element                              thon 12/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex27ScatraType Discret::ELEMENTS::SoHex27ScatraType::instance_;

Discret::ELEMENTS::SoHex27ScatraType& Discret::ELEMENTS::SoHex27ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoHex27ScatraType::Create(
    const std::vector<char>& data)
{
  auto* object =
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>(
          -1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex27ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex27ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoHex27ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex27;
  SoHex27Type::setup_element_definition(definitions_hex27);

  std::map<std::string, Input::LineDefinition>& defs_hex27 =
      definitions_hex27["SOLIDH27_DEPRECATED"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX27"] =
      Input::LineDefinition::Builder(defs_hex27["HEX27"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex27ScatraType::initialize(Core::FE::Discretization& dis)
{
  SoHex27Type::initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex27_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                       |
 *----------------------------------------------------------------------*/


Discret::ELEMENTS::SoTet4ScatraType Discret::ELEMENTS::SoTet4ScatraType::instance_;

Discret::ELEMENTS::SoTet4ScatraType& Discret::ELEMENTS::SoTet4ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoTet4ScatraType::Create(
    const std::vector<char>& data)
{
  auto* object =
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(-1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoTet4ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_tet4;
  SoTet4Type::setup_element_definition(definitions_tet4);

  std::map<std::string, Input::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4_DEPRECATED"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = Input::LineDefinition::Builder(defs_tet4["TET4"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoTet4ScatraType::initialize(Core::FE::Discretization& dis)
{
  SoTet4Type::initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 10 Element                                       |
 *----------------------------------------------------------------------*/


Discret::ELEMENTS::SoTet10ScatraType Discret::ELEMENTS::SoTet10ScatraType::instance_;

Discret::ELEMENTS::SoTet10ScatraType& Discret::ELEMENTS::SoTet10ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoTet10ScatraType::Create(
    const std::vector<char>& data)
{
  auto* object =
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>(
          -1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet10ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet10ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoTet10ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_tet10;
  SoTet10Type::setup_element_definition(definitions_tet10);

  std::map<std::string, Input::LineDefinition>& defs_tet10 =
      definitions_tet10["SOLIDT10_DEPRECATED"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET10"] =
      Input::LineDefinition::Builder(defs_tet10["TET10"]).add_named_string("TYPE").build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoTet10ScatraType::initialize(Core::FE::Discretization& dis)
{
  SoTet10Type::initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet10_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  WEDGE 6 Element                                       |
 *----------------------------------------------------------------------*/


Discret::ELEMENTS::SoWeg6ScatraType Discret::ELEMENTS::SoWeg6ScatraType::instance_;

Discret::ELEMENTS::SoWeg6ScatraType& Discret::ELEMENTS::SoWeg6ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoWeg6ScatraType::Create(
    const std::vector<char>& data)
{
  auto* object =
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoWeg6, Core::FE::CellType::wedge6>(
          -1, -1);
  object->unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoWeg6ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoWeg6, Core::FE::CellType::wedge6>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoWeg6ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoWeg6, Core::FE::CellType::wedge6>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoWeg6ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_weg6;
  SoWeg6Type::setup_element_definition(definitions_weg6);

  std::map<std::string, Input::LineDefinition>& defs_weg6 = definitions_weg6["SOLIDW6_DEPRECATED"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["WEDGE6"] =
      Input::LineDefinition::Builder(defs_weg6["WEDGE6"]).add_named_string("TYPE").build();
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoWeg6ScatraType::initialize(Core::FE::Discretization& dis)
{
  SoWeg6Type::initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoWeg6, Core::FE::CellType::wedge6>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_weg6_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
