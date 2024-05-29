/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1_poro_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_w1_poro.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                       |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallQuad4PoroType DRT::ELEMENTS::WallQuad4PoroType::instance_;

DRT::ELEMENTS::WallQuad4PoroType& DRT::ELEMENTS::WallQuad4PoroType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::WallQuad4PoroType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallQuad4PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ4PORO")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallQuad4PoroType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::WallQuad4PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLQ4PORO"];

  defs["QUAD4"] = INPUT::LineDefinition::Builder(defs_wall["QUAD4"])
                      .add_optional_named_double_vector("POROANISODIR1", 2)
                      .add_optional_named_double_vector("POROANISODIR2", 2)
                      .add_optional_named_double_vector("POROANISONODALCOEFFS1", 4)
                      .add_optional_named_double_vector("POROANISONODALCOEFFS2", 4)
                      .Build();
}

int DRT::ELEMENTS::WallQuad4PoroType::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad4>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::WallQuad9PoroType DRT::ELEMENTS::WallQuad9PoroType::instance_;

DRT::ELEMENTS::WallQuad9PoroType& DRT::ELEMENTS::WallQuad9PoroType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::WallQuad9PoroType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallQuad9PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ9PORO")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallQuad9PoroType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad9>(id, owner));
  return ele;
}

void DRT::ELEMENTS::WallQuad9PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLQ9PORO"];

  defs["QUAD9"] = INPUT::LineDefinition::Builder(defs_wall["QUAD9"])
                      .add_optional_named_double_vector("POROANISODIR1", 2)
                      .add_optional_named_double_vector("POROANISODIR2", 2)
                      .Build();
}

int DRT::ELEMENTS::WallQuad9PoroType::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::quad9>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  NURBS 4 Element                                       |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallNurbs4PoroType DRT::ELEMENTS::WallNurbs4PoroType::instance_;

DRT::ELEMENTS::WallNurbs4PoroType& DRT::ELEMENTS::WallNurbs4PoroType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallNurbs4PoroType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::nurbs4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallNurbs4PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLN4PORO")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::nurbs4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallNurbs4PoroType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::nurbs4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::WallNurbs4PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLN4PORO"];

  defs["NURBS4"] = INPUT::LineDefinition::Builder(defs_wall["NURBS4"])
                       .add_optional_named_double_vector("POROANISODIR1", 2)
                       .add_optional_named_double_vector("POROANISODIR2", 2)
                       .Build();
}

int DRT::ELEMENTS::WallNurbs4PoroType::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::nurbs4>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  NURBS 9 Element                                       |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallNurbs9PoroType DRT::ELEMENTS::WallNurbs9PoroType::instance_;

DRT::ELEMENTS::WallNurbs9PoroType& DRT::ELEMENTS::WallNurbs9PoroType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallNurbs9PoroType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::nurbs9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallNurbs9PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLN9PORO")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::nurbs9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallNurbs9PoroType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::nurbs9>(id, owner));
  return ele;
}

void DRT::ELEMENTS::WallNurbs9PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLN9PORO"];

  defs["NURBS9"] = INPUT::LineDefinition::Builder(defs_wall["NURBS9"])
                       .add_optional_named_double_vector("POROANISODIR1", 2)
                       .add_optional_named_double_vector("POROANISODIR2", 2)
                       .Build();
}

int DRT::ELEMENTS::WallNurbs9PoroType::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::nurbs9>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                       |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallTri3PoroType DRT::ELEMENTS::WallTri3PoroType::instance_;

DRT::ELEMENTS::WallTri3PoroType& DRT::ELEMENTS::WallTri3PoroType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::WallTri3PoroType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallTri3PoroType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLT3PORO")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallTri3PoroType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::tri3>(id, owner));
  return ele;
}

void DRT::ELEMENTS::WallTri3PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLT3PORO"];

  defs["TRI3"] = INPUT::LineDefinition::Builder(defs_wall["TRI3"])
                     .add_optional_named_double_vector("POROANISODIR1", 2)
                     .add_optional_named_double_vector("POROANISODIR2", 2)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS1", 3)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS2", 3)
                     .Build();
}

int DRT::ELEMENTS::WallTri3PoroType::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele =
        dynamic_cast<DRT::ELEMENTS::Wall1Poro<CORE::FE::CellType::tri3>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->InitElement();
  }
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
