/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element (p1/mixed approach).

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1_poro_p1_eletypes.hpp"

#include "4C_fluid_ele_nullspace.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_w1_poro_p1.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                                      |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallQuad4PoroP1Type DRT::ELEMENTS::WallQuad4PoroP1Type::instance_;

DRT::ELEMENTS::WallQuad4PoroP1Type& DRT::ELEMENTS::WallQuad4PoroP1Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallQuad4PoroP1Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallQuad4PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ4POROP1")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallQuad4PoroP1Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::WallQuad4PoroP1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wallporo;
  WallQuad4PoroType::setup_element_definition(definitions_wallporo);

  std::map<std::string, INPUT::LineDefinition>& defs_wallporo = definitions_wallporo["WALLQ4PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLQ4POROP1"];

  defs["QUAD4"] = defs_wallporo["QUAD4"];
}

void DRT::ELEMENTS::WallQuad4PoroP1Type::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 3;
  nv = 2;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::WallQuad4PoroP1Type::ComputeNullSpace(
    CORE::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

int DRT::ELEMENTS::WallQuad4PoroP1Type::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad4>* actele =
        dynamic_cast<DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad4>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_PoroP1* failed");
    actele->InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                                      |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallQuad9PoroP1Type DRT::ELEMENTS::WallQuad9PoroP1Type::instance_;

DRT::ELEMENTS::WallQuad9PoroP1Type& DRT::ELEMENTS::WallQuad9PoroP1Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallQuad9PoroP1Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallQuad9PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ9POROP1")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallQuad9PoroP1Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad9>(id, owner));
  return ele;
}

void DRT::ELEMENTS::WallQuad9PoroP1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wallporo;
  WallQuad9PoroType::setup_element_definition(definitions_wallporo);

  std::map<std::string, INPUT::LineDefinition>& defs_wallporo = definitions_wallporo["WALLQ9PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLQ9POROP1"];

  defs["QUAD9"] = defs_wallporo["QUAD9"];
}

void DRT::ELEMENTS::WallQuad9PoroP1Type::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 3;
  nv = 2;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::WallQuad9PoroP1Type::ComputeNullSpace(
    CORE::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

int DRT::ELEMENTS::WallQuad9PoroP1Type::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad9>* actele =
        dynamic_cast<DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad9>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_PoroP1* failed");
    actele->InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                                       |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallTri3PoroP1Type DRT::ELEMENTS::WallTri3PoroP1Type::instance_;

DRT::ELEMENTS::WallTri3PoroP1Type& DRT::ELEMENTS::WallTri3PoroP1Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallTri3PoroP1Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::tri3>* object =
      new DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallTri3PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLT3POROP1")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::WallTri3PoroP1Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::tri3>(id, owner));
  return ele;
}

void DRT::ELEMENTS::WallTri3PoroP1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wallporo;
  WallTri3PoroType::setup_element_definition(definitions_wallporo);

  std::map<std::string, INPUT::LineDefinition>& defs_wallporo = definitions_wallporo["WALLT3PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLT3POROP1"];

  defs["TRI3"] = defs_wallporo["TRI3"];
}

void DRT::ELEMENTS::WallTri3PoroP1Type::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 3;
  nv = 2;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::WallTri3PoroP1Type::ComputeNullSpace(
    CORE::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

int DRT::ELEMENTS::WallTri3PoroP1Type::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::tri3>* actele =
        dynamic_cast<DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::tri3>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_PoroP1* failed");
    actele->InitElement();
  }
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
