/*----------------------------------------------------------------------*/
/*!
 \file so3_poro_scatra_eletypes.cpp

 \brief element types of the 3D solid-poro element including scatra functionality

 \level 2

 \maintainer Christoph Schmidt
             schmidt@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 15251
*----------------------------------------------------------------------*/

#include "so3_poro_scatra.H"
#include "so3_poro_scatra_eletypes.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8PoroScatraType DRT::ELEMENTS::So_hex8PoroScatraType::instance_;

DRT::ELEMENTS::So_hex8PoroScatraType& DRT::ELEMENTS::So_hex8PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_hex8PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* object =
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex8;
  So_hex8PoroType::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8POROSCATRA"];

  defs["HEX8"] = defs_hex8["HEX8"];

  // add scalar transport ImplType
  defs["HEX8"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet4PoroScatraType DRT::ELEMENTS::So_tet4PoroScatraType::instance_;

DRT::ELEMENTS::So_tet4PoroScatraType& DRT::ELEMENTS::So_tet4PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_tet4PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* object =
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT4POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet4PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_tet4;
  So_tet4PoroType::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT4POROSCATRA"];

  defs["TET4"] = defs_tet4["TET4"];

  // add scalar transport ImplType
  defs["TET4"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  HEX 27 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_hex27PoroScatraType DRT::ELEMENTS::So_hex27PoroScatraType::instance_;

DRT::ELEMENTS::So_hex27PoroScatraType& DRT::ELEMENTS::So_hex27PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_hex27PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>* object =
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH27POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex27PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex27;
  So_hex27PoroType::SetupElementDefinition(definitions_hex27);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH27POROSCATRA"];

  defs["HEX27"] = defs_hex27["HEX27"];

  // add scalar transport ImplType
  defs["HEX27"].AddNamedString("TYPE");
}


/*----------------------------------------------------------------------*
 |  TET 10 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet10PoroScatraType DRT::ELEMENTS::So_tet10PoroScatraType::instance_;

DRT::ELEMENTS::So_tet10PoroScatraType& DRT::ELEMENTS::So_tet10PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_tet10PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>* object =
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT10POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet10PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_tet10;
  So_tet10PoroType::SetupElementDefinition(definitions_tet10);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT10POROSCATRA"];

  defs["TET10"] = defs_tet10["TET10"];

  // add scalar transport ImplType
  defs["TET10"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  NURBS 27 Element                                      schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_nurbs27PoroScatraType DRT::ELEMENTS::So_nurbs27PoroScatraType::instance_;

DRT::ELEMENTS::So_nurbs27PoroScatraType& DRT::ELEMENTS::So_nurbs27PoroScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_nurbs27PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>* object =
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_nurbs27PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SONURBS27POROSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_nurbs27PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Poro_Scatra<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>(
          id, owner));
  return ele;
}

void DRT::ELEMENTS::So_nurbs27PoroScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_nurbs27;
  So_nurbs27PoroType::SetupElementDefinition(definitions_nurbs27);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_nurbs27 =
      definitions_nurbs27["SONURBS27PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SONURBS27POROSCATRA"];

  defs["NURBS27"] = defs_nurbs27["NURBS27"];

  // add scalar transport ImplType
  defs["NURBS27"].AddNamedString("TYPE");
}
