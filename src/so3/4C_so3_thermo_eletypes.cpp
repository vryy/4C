// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_so3_thermo_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_thermo.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *  HEX8 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 08/12 |
 *----------------------------------------------------------------------*/
Discret::Elements::SoHex8ThermoType Discret::Elements::SoHex8ThermoType::instance_;


/*----------------------------------------------------------------------*
 | access an instance of thermo type                                    |
 *----------------------------------------------------------------------*/
Discret::Elements::SoHex8ThermoType& Discret::Elements::SoHex8ThermoType::instance()
{
  return instance_;
}


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::Elements::SoHex8ThermoType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::So3Thermo<Discret::Elements::SoHex8, Core::FE::CellType::hex8>* object =
      new Discret::Elements::So3Thermo<Discret::Elements::SoHex8, Core::FE::CellType::hex8>(-1, -1);
  object->unpack(buffer);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::Elements::SoHex8ThermoType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8THERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::make_rcp<
        Discret::Elements::So3Thermo<Discret::Elements::SoHex8, Core::FE::CellType::hex8>>(

        id, owner);
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::Elements::SoHex8ThermoType::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::make_rcp<
      Discret::Elements::So3Thermo<Discret::Elements::SoHex8, Core::FE::CellType::hex8>>(

      id, owner);
  return ele;

}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     dano 08/12 |
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex8ThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex8;
  SoHex8Type::setup_element_definition(definitions_hex8);

  std::map<std::string, Input::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8_DEPRECATED"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = defs_hex8["HEX8"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
int Discret::Elements::SoHex8ThermoType::initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.num_my_col_elements(); ++i)
  {
    if (dis.l_col_element(i)->element_type() != *this) continue;

    Discret::Elements::So3Thermo<Discret::Elements::SoHex8, Core::FE::CellType::hex8>* actele =
        dynamic_cast<
            Discret::Elements::So3Thermo<Discret::Elements::SoHex8, Core::FE::CellType::hex8>*>(
            dis.l_col_element(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8_thermo* failed");
    // initialise all quantities
    actele->SoHex8::init_jacobian_mapping();
    // as an alternative we can call: So_hex8Type::initialize(dis);
    actele->So3Thermo<Discret::Elements::SoHex8,
        Core::FE::CellType::hex8>::init_jacobian_mapping_special_for_tsi_elements(dis);
  }

  return 0;
}  // initialize()
/*----------------------------------------------------------------------------*
 | ENDE HEX8 Element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
