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

Discret::ELEMENTS::WallQuad4PoroP1Type Discret::ELEMENTS::WallQuad4PoroP1Type::instance_;

Discret::ELEMENTS::WallQuad4PoroP1Type& Discret::ELEMENTS::WallQuad4PoroP1Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallQuad4PoroP1Type::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad4PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ4POROP1")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad4PoroP1Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad4>(id, owner));
  return ele;
}

void Discret::ELEMENTS::WallQuad4PoroP1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wallporo;
  WallQuad4PoroType::setup_element_definition(definitions_wallporo);

  std::map<std::string, Input::LineDefinition>& defs_wallporo = definitions_wallporo["WALLQ4PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLQ4POROP1"];

  defs["QUAD4"] = defs_wallporo["QUAD4"];
}

void Discret::ELEMENTS::WallQuad4PoroP1Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 3;
  nv = 2;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::WallQuad4PoroP1Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

int Discret::ELEMENTS::WallQuad4PoroP1Type::Initialize(Discret::Discretization& dis)
{
  Discret::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad4>* actele =
        dynamic_cast<Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad4>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_PoroP1* failed");
    actele->InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                                      |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::WallQuad9PoroP1Type Discret::ELEMENTS::WallQuad9PoroP1Type::instance_;

Discret::ELEMENTS::WallQuad9PoroP1Type& Discret::ELEMENTS::WallQuad9PoroP1Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallQuad9PoroP1Type::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad9PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ9POROP1")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallQuad9PoroP1Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad9>(id, owner));
  return ele;
}

void Discret::ELEMENTS::WallQuad9PoroP1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wallporo;
  WallQuad9PoroType::setup_element_definition(definitions_wallporo);

  std::map<std::string, Input::LineDefinition>& defs_wallporo = definitions_wallporo["WALLQ9PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLQ9POROP1"];

  defs["QUAD9"] = defs_wallporo["QUAD9"];
}

void Discret::ELEMENTS::WallQuad9PoroP1Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 3;
  nv = 2;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::WallQuad9PoroP1Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

int Discret::ELEMENTS::WallQuad9PoroP1Type::Initialize(Discret::Discretization& dis)
{
  Discret::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad9>* actele =
        dynamic_cast<Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad9>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_PoroP1* failed");
    actele->InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                                       |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::WallTri3PoroP1Type Discret::ELEMENTS::WallTri3PoroP1Type::instance_;

Discret::ELEMENTS::WallTri3PoroP1Type& Discret::ELEMENTS::WallTri3PoroP1Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::WallTri3PoroP1Type::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::tri3>* object =
      new Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallTri3PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLT3POROP1")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::WallTri3PoroP1Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::tri3>(id, owner));
  return ele;
}

void Discret::ELEMENTS::WallTri3PoroP1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wallporo;
  WallTri3PoroType::setup_element_definition(definitions_wallporo);

  std::map<std::string, Input::LineDefinition>& defs_wallporo = definitions_wallporo["WALLT3PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLT3POROP1"];

  defs["TRI3"] = defs_wallporo["TRI3"];
}

void Discret::ELEMENTS::WallTri3PoroP1Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 3;
  nv = 2;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::WallTri3PoroP1Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

int Discret::ELEMENTS::WallTri3PoroP1Type::Initialize(Discret::Discretization& dis)
{
  Discret::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::tri3>* actele =
        dynamic_cast<Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::tri3>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_PoroP1* failed");
    actele->InitElement();
  }
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
