/*----------------------------------------------------------------------*/
/*! \file

\brief element types of the 3D solid-poro element


\level 2

*----------------------------------------------------------------------*/

#include "4C_so3_poro_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_poro.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8PoroType DRT::ELEMENTS::SoHex8PoroType::instance_;

DRT::ELEMENTS::SoHex8PoroType& DRT::ELEMENTS::SoHex8PoroType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex8PoroType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex8PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex8PoroType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoHex8PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8;
  SoHex8Type::setup_element_definition(definitions_hex8);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"])
                     .add_optional_named_double_vector("POROANISODIR1", 3)
                     .add_optional_named_double_vector("POROANISODIR2", 3)
                     .add_optional_named_double_vector("POROANISODIR3", 3)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS1", 8)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS2", 8)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS3", 8)
                     .Build();
}

int DRT::ELEMENTS::SoHex8PoroType::Initialize(DRT::Discretization& dis)
{
  SoHex8Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8_poro* failed");
    actele->InitElement();
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  TET 4 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoTet4PoroType DRT::ELEMENTS::SoTet4PoroType::instance_;

DRT::ELEMENTS::SoTet4PoroType& DRT::ELEMENTS::SoTet4PoroType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoTet4PoroType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoTet4PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoTet4PoroType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoTet4PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet4;
  SoTet4Type::setup_element_definition(definitions_tet4);

  std::map<std::string, INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = INPUT::LineDefinition::Builder(defs_tet4["TET4"])
                     .add_optional_named_double_vector("POROANISODIR1", 3)
                     .add_optional_named_double_vector("POROANISODIR2", 3)
                     .add_optional_named_double_vector("POROANISODIR3", 3)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS1", 4)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS2", 4)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS3", 4)
                     .Build();
}

int DRT::ELEMENTS::SoTet4PoroType::Initialize(DRT::Discretization& dis)
{
  SoTet4Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4_poro* failed");
    actele->So3Poro<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>::InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  HEX 27 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex27PoroType DRT::ELEMENTS::SoHex27PoroType::instance_;

DRT::ELEMENTS::SoHex27PoroType& DRT::ELEMENTS::SoHex27PoroType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex27PoroType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex27PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex27PoroType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoHex27PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex27;
  SoHex27Type::setup_element_definition(definitions_hex27);

  std::map<std::string, INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX27"] = INPUT::LineDefinition::Builder(defs_hex27["HEX27"])
                      .add_optional_named_double_vector("POROANISODIR1", 3)
                      .add_optional_named_double_vector("POROANISODIR2", 3)
                      .add_optional_named_double_vector("POROANISODIR3", 3)
                      .Build();
}

int DRT::ELEMENTS::SoHex27PoroType::Initialize(DRT::Discretization& dis)
{
  SoHex27Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex27_poro* failed");
    actele->So3Poro<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>::InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 10 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoTet10PoroType DRT::ELEMENTS::SoTet10PoroType::instance_;

DRT::ELEMENTS::SoTet10PoroType& DRT::ELEMENTS::SoTet10PoroType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoTet10PoroType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoTet10PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoTet10PoroType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoTet10PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet10;
  SoTet10Type::setup_element_definition(definitions_tet10);

  std::map<std::string, INPUT::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET10"] = INPUT::LineDefinition::Builder(defs_tet10["TET10"])
                      .add_optional_named_double_vector("POROANISODIR1", 3)
                      .add_optional_named_double_vector("POROANISODIR2", 3)
                      .add_optional_named_double_vector("POROANISODIR3", 3)
                      .Build();
}

int DRT::ELEMENTS::SoTet10PoroType::Initialize(DRT::Discretization& dis)
{
  SoTet10Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet10_poro* failed");
    actele->So3Poro<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>::InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  NURBS 27 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoNurbs27PoroType DRT::ELEMENTS::SoNurbs27PoroType::instance_;

DRT::ELEMENTS::SoNurbs27PoroType& DRT::ELEMENTS::SoNurbs27PoroType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::SoNurbs27PoroType::Create(const std::vector<char>& data)
{
  auto* object =
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::NURBS::SoNurbs27, CORE::FE::CellType::nurbs27>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoNurbs27PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::NURBS::SoNurbs27, CORE::FE::CellType::nurbs27>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoNurbs27PoroType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::NURBS::SoNurbs27, CORE::FE::CellType::nurbs27>(
          id, owner));
  return ele;
}

void DRT::ELEMENTS::SoNurbs27PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_nurbs27;
  NURBS::SoNurbs27Type::setup_element_definition(definitions_nurbs27);

  std::map<std::string, INPUT::LineDefinition>& defs_nurbs27 = definitions_nurbs27["SONURBS27"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["NURBS27"] = INPUT::LineDefinition::Builder(defs_nurbs27["NURBS27"])
                        .add_optional_named_double_vector("POROANISODIR1", 3)
                        .add_optional_named_double_vector("POROANISODIR2", 3)
                        .add_optional_named_double_vector("POROANISODIR3", 3)
                        .Build();
}

int DRT::ELEMENTS::SoNurbs27PoroType::Initialize(DRT::Discretization& dis)
{
  NURBS::SoNurbs27Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        DRT::ELEMENTS::So3Poro<DRT::ELEMENTS::NURBS::SoNurbs27, CORE::FE::CellType::nurbs27>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_nurbs27_poro* failed");
    actele->So3Poro<DRT::ELEMENTS::NURBS::SoNurbs27, CORE::FE::CellType::nurbs27>::InitElement();
  }
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
