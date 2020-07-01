/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D wall element for solid-part of porous medium using p1 (mixed)
 approach including scatra functionality.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "wall1_poro_p1_scatra_eletypes.H"
#include "wall1_poro_p1_scatra.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallQuad4PoroP1ScatraType DRT::ELEMENTS::WallQuad4PoroP1ScatraType::instance_;

DRT::ELEMENTS::WallQuad4PoroP1ScatraType& DRT::ELEMENTS::WallQuad4PoroP1ScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::WallQuad4PoroP1ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad4>* object =
      new DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad4PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ4POROP1SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad4PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad4>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallQuad4PoroP1ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_wallporo;
  WallQuad4PoroP1Type::SetupElementDefinition(definitions_wallporo);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_wallporo =
      definitions_wallporo["WALLQ4POROP1"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["WALLQ4POROP1SCATRA"];

  defs["QUAD4"] = defs_wallporo["QUAD4"];

  // add scalar transport ImplType
  defs["QUAD4"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallQuad9PoroP1ScatraType DRT::ELEMENTS::WallQuad9PoroP1ScatraType::instance_;

DRT::ELEMENTS::WallQuad9PoroP1ScatraType& DRT::ELEMENTS::WallQuad9PoroP1ScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::WallQuad9PoroP1ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad9>* object =
      new DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad9PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ9POROP1SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad9PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad9>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallQuad9PoroP1ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_wallporo;
  WallQuad9PoroP1Type::SetupElementDefinition(definitions_wallporo);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_wallporo =
      definitions_wallporo["WALLQ9POROP1"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["WALLQ9POROP1SCATRA"];

  defs["QUAD9"] = defs_wallporo["QUAD9"];

  // add scalar transport ImplType
  defs["QUAD9"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallTri3PoroP1ScatraType DRT::ELEMENTS::WallTri3PoroP1ScatraType::instance_;

DRT::ELEMENTS::WallTri3PoroP1ScatraType& DRT::ELEMENTS::WallTri3PoroP1ScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::WallTri3PoroP1ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::tri3>* object =
      new DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallTri3PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLT3POROP1SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallTri3PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::tri3>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallTri3PoroP1ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_wallporo;
  WallTri3PoroP1Type::SetupElementDefinition(definitions_wallporo);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_wallporo =
      definitions_wallporo["WALLT3POROP1"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["WALLT3POROP1SCATRA"];

  defs["TRI3"] = defs_wallporo["TRI3"];

  // add scalar transport ImplType
  defs["TRI3"].AddNamedString("TYPE");
}
