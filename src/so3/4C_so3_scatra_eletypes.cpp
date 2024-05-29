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
DRT::ELEMENTS::SoHex8ScatraType DRT::ELEMENTS::SoHex8ScatraType::instance_;

DRT::ELEMENTS::SoHex8ScatraType& DRT::ELEMENTS::SoHex8ScatraType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex8ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex8ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex8ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoHex8ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8;
  SoHex8Type::setup_element_definition(definitions_hex8);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex8ScatraType::Initialize(DRT::Discretization& dis)
{
  SoHex8Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8_scatra* failed");
    actele->InitElement();
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  HEX 8 fbar Element                                       Thon 12/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8fbarScatraType DRT::ELEMENTS::SoHex8fbarScatraType::instance_;

DRT::ELEMENTS::SoHex8fbarScatraType& DRT::ELEMENTS::SoHex8fbarScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex8fbarScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8fbar, CORE::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex8fbarScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8fbar, CORE::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex8fbarScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8fbar, CORE::FE::CellType::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoHex8fbarScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8;
  SoHex8fbarType::setup_element_definition(definitions_hex8);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8FBAR"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex8fbarScatraType::Initialize(DRT::Discretization& dis)
{
  SoHex8fbarType::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex8fbar, CORE::FE::CellType::hex8>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8fbar_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  HEX 27 Solid Scatra Element                              thon 12/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex27ScatraType DRT::ELEMENTS::SoHex27ScatraType::instance_;

DRT::ELEMENTS::SoHex27ScatraType& DRT::ELEMENTS::SoHex27ScatraType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex27ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex27ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex27ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoHex27ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex27;
  SoHex27Type::setup_element_definition(definitions_hex27);

  std::map<std::string, INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX27"] =
      INPUT::LineDefinition::Builder(defs_hex27["HEX27"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex27ScatraType::Initialize(DRT::Discretization& dis)
{
  SoHex27Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex27_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::SoTet4ScatraType DRT::ELEMENTS::SoTet4ScatraType::instance_;

DRT::ELEMENTS::SoTet4ScatraType& DRT::ELEMENTS::SoTet4ScatraType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoTet4ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet4ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet4ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoTet4ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet4;
  SoTet4Type::setup_element_definition(definitions_tet4);

  std::map<std::string, INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = INPUT::LineDefinition::Builder(defs_tet4["TET4"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoTet4ScatraType::Initialize(DRT::Discretization& dis)
{
  SoTet4Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 10 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::SoTet10ScatraType DRT::ELEMENTS::SoTet10ScatraType::instance_;

DRT::ELEMENTS::SoTet10ScatraType& DRT::ELEMENTS::SoTet10ScatraType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoTet10ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet10ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet10ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoTet10ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet10;
  SoTet10Type::setup_element_definition(definitions_tet10);

  std::map<std::string, INPUT::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET10"] =
      INPUT::LineDefinition::Builder(defs_tet10["TET10"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoTet10ScatraType::Initialize(DRT::Discretization& dis)
{
  SoTet10Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet10_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  WEDGE 6 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::SoWeg6ScatraType DRT::ELEMENTS::SoWeg6ScatraType::instance_;

DRT::ELEMENTS::SoWeg6ScatraType& DRT::ELEMENTS::SoWeg6ScatraType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoWeg6ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoWeg6, CORE::FE::CellType::wedge6>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoWeg6ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoWeg6, CORE::FE::CellType::wedge6>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoWeg6ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoWeg6, CORE::FE::CellType::wedge6>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoWeg6ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_weg6;
  SoWeg6Type::setup_element_definition(definitions_weg6);

  std::map<std::string, INPUT::LineDefinition>& defs_weg6 = definitions_weg6["SOLIDW6"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["WEDGE6"] =
      INPUT::LineDefinition::Builder(defs_weg6["WEDGE6"]).AddNamedString("TYPE").Build();
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoWeg6ScatraType::Initialize(DRT::Discretization& dis)
{
  SoWeg6Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Scatra<DRT::ELEMENTS::SoWeg6, CORE::FE::CellType::wedge6>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_weg6_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
