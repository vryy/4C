/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element including scatra functionality.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "baci_w1_poro_scatra_eletypes.H"

#include "baci_io_linedefinition.H"
#include "baci_w1_poro_scatra.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallQuad4PoroScatraType DRT::ELEMENTS::WallQuad4PoroScatraType::instance_;

DRT::ELEMENTS::WallQuad4PoroScatraType& DRT::ELEMENTS::WallQuad4PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallQuad4PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad4>* object =
      new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad4PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ4POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad4PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad4>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallQuad4PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  WallQuad4PoroType::SetupElementDefinition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALLQ4PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLQ4POROSCATRA"];

  defs["QUAD4"] = INPUT::LineDefinition::Builder(defs_wall["QUAD4"]).AddNamedString("TYPE").Build();
}


/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::WallQuad9PoroScatraType DRT::ELEMENTS::WallQuad9PoroScatraType::instance_;

DRT::ELEMENTS::WallQuad9PoroScatraType& DRT::ELEMENTS::WallQuad9PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallQuad9PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad9>* object =
      new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad9PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ9POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad9PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::quad9>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallQuad9PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  WallQuad9PoroType::SetupElementDefinition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALLQ9PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLQ9POROSCATRA"];

  defs["QUAD9"] = INPUT::LineDefinition::Builder(defs_wall["QUAD9"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  NURBS 4 Element                                       schmidt 09/17 |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallNurbs4PoroScatraType DRT::ELEMENTS::WallNurbs4PoroScatraType::instance_;

DRT::ELEMENTS::WallNurbs4PoroScatraType& DRT::ELEMENTS::WallNurbs4PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallNurbs4PoroScatraType::Create(
    const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs4>* object =
      new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs4>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallNurbs4PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLN4POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallNurbs4PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs4>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallNurbs4PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  WallNurbs4PoroType::SetupElementDefinition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALLN4PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLN4POROSCATRA"];

  defs["NURBS4"] =
      INPUT::LineDefinition::Builder(defs_wall["NURBS4"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  NURBS 9 Element                                       schmidt 09/17 |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallNurbs9PoroScatraType DRT::ELEMENTS::WallNurbs9PoroScatraType::instance_;

DRT::ELEMENTS::WallNurbs9PoroScatraType& DRT::ELEMENTS::WallNurbs9PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallNurbs9PoroScatraType::Create(
    const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs9>* object =
      new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs9>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallNurbs9PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLN9POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallNurbs9PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::nurbs9>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallNurbs9PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  WallNurbs9PoroType::SetupElementDefinition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALLN9PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLN9POROSCATRA"];

  defs["NURBS9"] =
      INPUT::LineDefinition::Builder(defs_wall["NURBS9"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallTri3PoroScatraType DRT::ELEMENTS::WallTri3PoroScatraType::instance_;

DRT::ELEMENTS::WallTri3PoroScatraType& DRT::ELEMENTS::WallTri3PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::WallTri3PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::tri3>* object =
      new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallTri3PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLT3POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallTri3PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Wall1_Poro_Scatra<CORE::FE::CellType::tri3>(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallTri3PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_wall;
  WallTri3PoroType::SetupElementDefinition(definitions_wall);

  std::map<std::string, INPUT::LineDefinition>& defs_wall = definitions_wall["WALLT3PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["WALLT3POROSCATRA"];

  defs["TRI3"] = INPUT::LineDefinition::Builder(defs_wall["TRI3"]).AddNamedString("TYPE").Build();
}

BACI_NAMESPACE_CLOSE
