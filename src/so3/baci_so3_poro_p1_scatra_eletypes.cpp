/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the 3D solid-poro element (p1, mixed approach) including scatra
 functionality

 \level 2

 *----------------------------------------------------------------------*/

#include "baci_so3_poro_p1_scatra_eletypes.hpp"

#include "baci_io_linedefinition.hpp"
#include "baci_so3_poro_p1_scatra.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8PoroP1ScatraType DRT::ELEMENTS::So_hex8PoroP1ScatraType::instance_;

DRT::ELEMENTS::So_hex8PoroP1ScatraType& DRT::ELEMENTS::So_hex8PoroP1ScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::So_hex8PoroP1ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>* object =
      new DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_hex8, CORE::FE::CellType::hex8>(
          id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8PoroP1ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8poro;
  So_hex8PoroP1Type::SetupElementDefinition(definitions_hex8poro);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8poro["SOLIDH8POROP1"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet4PoroP1ScatraType DRT::ELEMENTS::So_tet4PoroP1ScatraType::instance_;

DRT::ELEMENTS::So_tet4PoroP1ScatraType& DRT::ELEMENTS::So_tet4PoroP1ScatraType::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::So_tet4PoroP1ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>* object =
      new DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_tet4, CORE::FE::CellType::tet4>(
          id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet4PoroP1ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet4;
  So_tet4PoroP1Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4POROP1"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["TET4"] = INPUT::LineDefinition::Builder(defs_tet4["TET4"]).AddNamedString("TYPE").Build();
}

BACI_NAMESPACE_CLOSE
