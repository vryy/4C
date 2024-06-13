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
Discret::ELEMENTS::SoHex8PoroP1ScatraType Discret::ELEMENTS::SoHex8PoroP1ScatraType::instance_;

Discret::ELEMENTS::SoHex8PoroP1ScatraType& Discret::ELEMENTS::SoHex8PoroP1ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoHex8PoroP1ScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>* object =
      new Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoHex8PoroP1ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex8poro;
  SoHex8PoroP1Type::setup_element_definition(definitions_hex8poro);

  std::map<std::string, Input::LineDefinition>& defs_hex8 = definitions_hex8poro["SOLIDH8POROP1"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = Input::LineDefinition::Builder(defs_hex8["HEX8"]).add_named_string("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoTet4PoroP1ScatraType Discret::ELEMENTS::SoTet4PoroP1ScatraType::instance_;

Discret::ELEMENTS::SoTet4PoroP1ScatraType& Discret::ELEMENTS::SoTet4PoroP1ScatraType::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoTet4PoroP1ScatraType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>* object =
      new Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
          -1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4PoroP1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4PoroP1ScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoTet4PoroP1ScatraType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_tet4;
  SoTet4PoroP1Type::setup_element_definition(definitions_tet4);

  std::map<std::string, Input::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4POROP1"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = Input::LineDefinition::Builder(defs_tet4["TET4"]).add_named_string("TYPE").Build();
}

FOUR_C_NAMESPACE_CLOSE
