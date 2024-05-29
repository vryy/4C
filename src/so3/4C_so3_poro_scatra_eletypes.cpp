/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the 3D solid-poro element including scatra functionality

 \level 2

*----------------------------------------------------------------------*/

#include "4C_so3_poro_scatra_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_poro_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8PoroScatraType DRT::ELEMENTS::SoHex8PoroScatraType::instance_;

DRT::ELEMENTS::SoHex8PoroScatraType& DRT::ELEMENTS::SoHex8PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex8PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>* object =
      new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex8PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex8PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoHex8PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8;
  SoHex8PoroType::setup_element_definition(definitions_hex8);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::SoTet4PoroScatraType DRT::ELEMENTS::SoTet4PoroScatraType::instance_;

DRT::ELEMENTS::SoTet4PoroScatraType& DRT::ELEMENTS::SoTet4PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::SoTet4PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>* object =
      new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet4PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet4PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoTet4PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet4;
  SoTet4PoroType::setup_element_definition(definitions_tet4);

  std::map<std::string, INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = INPUT::LineDefinition::Builder(defs_tet4["TET4"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  HEX 27 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::SoHex27PoroScatraType DRT::ELEMENTS::SoHex27PoroScatraType::instance_;

DRT::ELEMENTS::SoHex27PoroScatraType& DRT::ELEMENTS::SoHex27PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex27PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>* object =
      new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex27PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex27PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>(
          id, owner));
  return ele;
}

void DRT::ELEMENTS::SoHex27PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex27;
  SoHex27PoroType::setup_element_definition(definitions_hex27);

  std::map<std::string, INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX27"] =
      INPUT::LineDefinition::Builder(defs_hex27["HEX27"]).AddNamedString("TYPE").Build();
}


/*----------------------------------------------------------------------*
 |  TET 10 Element                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::SoTet10PoroScatraType DRT::ELEMENTS::SoTet10PoroScatraType::instance_;

DRT::ELEMENTS::SoTet10PoroScatraType& DRT::ELEMENTS::SoTet10PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::SoTet10PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>* object =
      new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet10PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet10PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>(
          id, owner));
  return ele;
}

void DRT::ELEMENTS::SoTet10PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet10;
  SoTet10PoroType::setup_element_definition(definitions_tet10);

  std::map<std::string, INPUT::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET10"] =
      INPUT::LineDefinition::Builder(defs_tet10["TET10"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  NURBS 27 Element                                      schmidt 09/17 |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::SoNurbs27PoroScatraType DRT::ELEMENTS::SoNurbs27PoroScatraType::instance_;

DRT::ELEMENTS::SoNurbs27PoroScatraType& DRT::ELEMENTS::SoNurbs27PoroScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::SoNurbs27PoroScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::NURBS::SoNurbs27, CORE::FE::CellType::nurbs27>*
      object = new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::NURBS::SoNurbs27,
          CORE::FE::CellType::nurbs27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoNurbs27PoroScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::NURBS::SoNurbs27,
            CORE::FE::CellType::nurbs27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoNurbs27PoroScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::NURBS::SoNurbs27,
          CORE::FE::CellType::nurbs27>(id, owner));
  return ele;
}

void DRT::ELEMENTS::SoNurbs27PoroScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_nurbs27;
  SoNurbs27PoroType::setup_element_definition(definitions_nurbs27);

  std::map<std::string, INPUT::LineDefinition>& defs_nurbs27 = definitions_nurbs27["SONURBS27PORO"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["NURBS27"] =
      INPUT::LineDefinition::Builder(defs_nurbs27["NURBS27"]).AddNamedString("TYPE").Build();
}

FOUR_C_NAMESPACE_CLOSE
