/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the 3D solid-poro element (p1, mixed approach) including scatra
 functionality

 \level 2

 *----------------------------------------------------------------------*/

#include "4C_so3_poro_p1_scatra_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_poro_p1_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8PoroP1ScatraType DRT::ELEMENTS::SoHex8PoroP1ScatraType::instance_;

DRT::ELEMENTS::SoHex8PoroP1ScatraType& DRT::ELEMENTS::SoHex8PoroP1ScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::SoHex8PoroP1ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>* object =
      new DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex8PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoHex8PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>(
          id, owner));
  return ele;
}

void DRT::ELEMENTS::SoHex8PoroP1ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8poro;
  SoHex8PoroP1Type::setup_element_definition(definitions_hex8poro);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8poro["SOLIDH8POROP1"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoTet4PoroP1ScatraType DRT::ELEMENTS::SoTet4PoroP1ScatraType::instance_;

DRT::ELEMENTS::SoTet4PoroP1ScatraType& DRT::ELEMENTS::SoTet4PoroP1ScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::SoTet4PoroP1ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>* object =
      new DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet4PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet4PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3PoroP1Scatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>(
          id, owner));
  return ele;
}

void DRT::ELEMENTS::SoTet4PoroP1ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet4;
  SoTet4PoroP1Type::setup_element_definition(definitions_tet4);

  std::map<std::string, INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4POROP1"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = INPUT::LineDefinition::Builder(defs_tet4["TET4"]).AddNamedString("TYPE").Build();
}

FOUR_C_NAMESPACE_CLOSE
