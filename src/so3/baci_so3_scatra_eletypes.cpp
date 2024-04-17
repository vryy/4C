/*----------------------------------------------------------------------*/
/*! \file

\brief Solid-scatra elements types

\level 2


*----------------------------------------------------------------------*/

#include "baci_io_linedefinition.hpp"
#include "baci_so3_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8ScatraType DRT::ELEMENTS::So_hex8ScatraType::instance_;

DRT::ELEMENTS::So_hex8ScatraType& DRT::ELEMENTS::So_hex8ScatraType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::So_hex8ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8;
  So_hex8Type::SetupElementDefinition(definitions_hex8);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8ScatraType::Initialize(DRT::Discretization& dis)
{
  So_hex8Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8_scatra* failed");
    actele->InitElement();
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  HEX 8 fbar Element                                       Thon 12/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8fbarScatraType DRT::ELEMENTS::So_hex8fbarScatraType::instance_;

DRT::ELEMENTS::So_hex8fbarScatraType& DRT::ELEMENTS::So_hex8fbarScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::So_hex8fbarScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, CORE::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, CORE::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, CORE::FE::CellType::hex8>(
          id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8fbarScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8;
  So_hex8fbarType::SetupElementDefinition(definitions_hex8);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8FBAR"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbarScatraType::Initialize(DRT::Discretization& dis)
{
  So_hex8fbarType::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, CORE::FE::CellType::hex8>*>(
        dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8fbar_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  HEX 27 Solid Scatra Element                              thon 12/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex27ScatraType DRT::ELEMENTS::So_hex27ScatraType::instance_;

DRT::ELEMENTS::So_hex27ScatraType& DRT::ELEMENTS::So_hex27ScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::So_hex27ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, CORE::FE::CellType::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, CORE::FE::CellType::hex27>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, CORE::FE::CellType::hex27>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex27ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex27;
  So_hex27Type::SetupElementDefinition(definitions_hex27);

  std::map<std::string, INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX27"] =
      INPUT::LineDefinition::Builder(defs_hex27["HEX27"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex27ScatraType::Initialize(DRT::Discretization& dis)
{
  So_hex27Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, CORE::FE::CellType::hex27>*>(
        dis.lColElement(i));
    if (!actele) dserror("cast to So_hex27_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet4ScatraType DRT::ELEMENTS::So_tet4ScatraType::instance_;

DRT::ELEMENTS::So_tet4ScatraType& DRT::ELEMENTS::So_tet4ScatraType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::So_tet4ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet4ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet4;
  So_tet4Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["TET4"] = INPUT::LineDefinition::Builder(defs_tet4["TET4"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4ScatraType::Initialize(DRT::Discretization& dis)
{
  So_tet4Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 10 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet10ScatraType DRT::ELEMENTS::So_tet10ScatraType::instance_;

DRT::ELEMENTS::So_tet10ScatraType& DRT::ELEMENTS::So_tet10ScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::So_tet10ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, CORE::FE::CellType::tet10>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, CORE::FE::CellType::tet10>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, CORE::FE::CellType::tet10>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet10ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet10;
  So_tet10Type::SetupElementDefinition(definitions_tet10);

  std::map<std::string, INPUT::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["TET10"] =
      INPUT::LineDefinition::Builder(defs_tet10["TET10"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10ScatraType::Initialize(DRT::Discretization& dis)
{
  So_tet10Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, CORE::FE::CellType::tet10>*>(
        dis.lColElement(i));
    if (!actele) dserror("cast to So_tet10_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  WEDGE 6 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_weg6ScatraType DRT::ELEMENTS::So_weg6ScatraType::instance_;

DRT::ELEMENTS::So_weg6ScatraType& DRT::ELEMENTS::So_weg6ScatraType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::So_weg6ScatraType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, CORE::FE::CellType::wedge6>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_weg6ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, CORE::FE::CellType::wedge6>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_weg6ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, CORE::FE::CellType::wedge6>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_weg6ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_weg6;
  So_weg6Type::SetupElementDefinition(definitions_weg6);

  std::map<std::string, INPUT::LineDefinition>& defs_weg6 = definitions_weg6["SOLIDW6"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["WEDGE6"] =
      INPUT::LineDefinition::Builder(defs_weg6["WEDGE6"]).AddNamedString("TYPE").Build();
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_weg6ScatraType::Initialize(DRT::Discretization& dis)
{
  So_weg6Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, CORE::FE::CellType::wedge6>*>(
        dis.lColElement(i));
    if (!actele) dserror("cast to So_weg6_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
