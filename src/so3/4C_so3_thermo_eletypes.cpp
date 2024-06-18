/*----------------------------------------------------------------------*/
/*! \file
\brief 3d TSI solid element
\level 1

*----------------------------------------------------------------------*/

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
Discret::ELEMENTS::SoHex8ThermoType Discret::ELEMENTS::SoHex8ThermoType::instance_;


/*----------------------------------------------------------------------*
 | access an instance of thermo type                                    |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8ThermoType& Discret::ELEMENTS::SoHex8ThermoType::Instance()
{
  return instance_;
}


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoHex8ThermoType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>* object =
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(-1, -1);
  object->unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8THERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8ThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
          id, owner));
  return ele;

}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     dano 08/12 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8ThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex8;
  SoHex8Type::setup_element_definition(definitions_hex8);

  std::map<std::string, Input::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = defs_hex8["HEX8"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex8ThermoType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>* actele =
        dynamic_cast<
            Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8_thermo* failed");
    // initialise all quantities
    actele->SoHex8::init_jacobian_mapping();
    // as an alternative we can call: So_hex8Type::Initialize(dis);
    actele->So3Thermo<Discret::ELEMENTS::SoHex8,
        Core::FE::CellType::hex8>::init_jacobian_mapping_special_for_nurbs(dis);
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE HEX8 Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *  HEX8FBAR element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 05/13 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8fbarThermoType Discret::ELEMENTS::SoHex8fbarThermoType::instance_;

Discret::ELEMENTS::SoHex8fbarThermoType& Discret::ELEMENTS::SoHex8fbarThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 05/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoHex8fbarThermoType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>* object =
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>(
          -1, -1);
  object->unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 05/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8fbarThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8FBARTHERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 05/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8fbarThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>(
          id, owner));
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     dano 05/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8fbarThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  // original definition
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex8fbar;

  // call setup of so3_ele
  SoHex8fbarType::setup_element_definition(definitions_hex8fbar);

  std::map<std::string, Input::LineDefinition>& defs_hex8fbar = definitions_hex8fbar["SOLIDH8FBAR"];

  // templated definition
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = defs_hex8fbar["HEX8"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex8fbarThermoType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>* actele =
        dynamic_cast<
            Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex8fbar, Core::FE::CellType::hex8>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8fbar_thermo* failed");

    // initialise all quantities
    actele->SoHex8fbar::init_jacobian_mapping();
    // as an alternative we can call: So_hex8fbarType::Initialize(dis);
    actele->So3Thermo<Discret::ELEMENTS::SoHex8fbar,
        Core::FE::CellType::hex8>::init_jacobian_mapping_special_for_nurbs(dis);
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE HEX8FBAR Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *  TET4 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 08/12 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoTet4ThermoType Discret::ELEMENTS::SoTet4ThermoType::instance_;

Discret::ELEMENTS::SoTet4ThermoType& Discret::ELEMENTS::SoTet4ThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoTet4ThermoType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>* object =
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(-1, -1);
  object->unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT4THERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4ThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
          id, owner));
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 08/12 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoTet4ThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_tet4;
  SoTet4Type::setup_element_definition(definitions_tet4);

  std::map<std::string, Input::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = defs_tet4["TET4"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoTet4ThermoType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>* actele =
        dynamic_cast<
            Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4_thermo* failed");

    actele->SoTet4::init_jacobian_mapping();
    // as an alternative we can call: So_tet4Type::Initialize(dis);
    actele->So3Thermo<Discret::ELEMENTS::SoTet4,
        Core::FE::CellType::tet4>::init_jacobian_mapping_special_for_nurbs(dis);
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE TET4 Element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 *  TET10 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                         farah 05/14 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoTet10ThermoType Discret::ELEMENTS::SoTet10ThermoType::instance_;

Discret::ELEMENTS::SoTet10ThermoType& Discret::ELEMENTS::SoTet10ThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoTet10ThermoType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>* object =
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>(
          -1, -1);
  object->unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet10ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT10THERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet10ThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>(
          id, owner));
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | build an instance of thermo type                         farah 05/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoTet10ThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_tet10;
  SoTet10Type::setup_element_definition(definitions_tet10);

  std::map<std::string, Input::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET10"] = defs_tet10["TET10"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          farah 05/14 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoTet10ThermoType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>* actele =
        dynamic_cast<
            Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet10_thermo* failed");

    actele->SoTet10::init_jacobian_mapping();
    // as an alternative we can call: So_tet4Type::Initialize(dis);
    actele->So3Thermo<Discret::ELEMENTS::SoTet10,
        Core::FE::CellType::tet10>::init_jacobian_mapping_special_for_nurbs(dis);
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE TET10 Element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  HEX 27 Element
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 10/13 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex27ThermoType Discret::ELEMENTS::SoHex27ThermoType::instance_;

Discret::ELEMENTS::SoHex27ThermoType& Discret::ELEMENTS::SoHex27ThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 10/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoHex27ThermoType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>* object =
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>(
          -1, -1);
  object->unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 10/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex27ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH27THERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 10/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex27ThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>(
          id, owner));
  return ele;
}  // Create ()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     dano 10/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex27ThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex27;
  SoHex27Type::setup_element_definition(definitions_hex27);

  std::map<std::string, Input::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX27"] = defs_hex27["HEX27"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 10/13 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex27ThermoType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>* actele =
        dynamic_cast<
            Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex27_thermo* failed");

    actele->SoHex27::init_jacobian_mapping();
    // as an alternative we can call: So_hex27Type::Initialize(dis);
    actele->So3Thermo<Discret::ELEMENTS::SoHex27,
        Core::FE::CellType::hex27>::init_jacobian_mapping_special_for_nurbs(dis);
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE HEX27 Element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  HEX 20 Element
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                         farah 05/14 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex20ThermoType Discret::ELEMENTS::SoHex20ThermoType::instance_;

Discret::ELEMENTS::SoHex20ThermoType& Discret::ELEMENTS::SoHex20ThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoHex20ThermoType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex20, Core::FE::CellType::hex20>* object =
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex20, Core::FE::CellType::hex20>(
          -1, -1);
  object->unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex20ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH20THERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex20, Core::FE::CellType::hex20>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex20ThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex20, Core::FE::CellType::hex20>(
          id, owner));
  return ele;
}  // Create ()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                    farah 05/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex20ThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex20;
  SoHex20Type::setup_element_definition(definitions_hex20);

  std::map<std::string, Input::LineDefinition>& defs_hex20 = definitions_hex20["SOLIDH20"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX20"] = defs_hex20["HEX20"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          farah 05/14 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex20ThermoType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex20, Core::FE::CellType::hex20>* actele =
        dynamic_cast<
            Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::SoHex20, Core::FE::CellType::hex20>*>(
            dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex20_thermo* failed");

    actele->SoHex20::init_jacobian_mapping();
    // as an alternative we can call: So_hex27Type::Initialize(dis);
    actele->So3Thermo<Discret::ELEMENTS::SoHex20,
        Core::FE::CellType::hex20>::init_jacobian_mapping_special_for_nurbs(dis);
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | ENDE HEX20 Element
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 |  nurbs 27 Element
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                         seitz 12/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoNurbs27ThermoType Discret::ELEMENTS::SoNurbs27ThermoType::instance_;

Discret::ELEMENTS::SoNurbs27ThermoType& Discret::ELEMENTS::SoNurbs27ThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 12/15 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoNurbs27ThermoType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::Nurbs::SoNurbs27, Core::FE::CellType::nurbs27>*
      object = new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::Nurbs::SoNurbs27,
          Core::FE::CellType::nurbs27>(-1, -1);
  object->unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 12/15 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoNurbs27ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SONURBS27THERMO")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::Nurbs::SoNurbs27,
            Core::FE::CellType::nurbs27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 12/15 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoNurbs27ThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::Nurbs::SoNurbs27,
          Core::FE::CellType::nurbs27>(id, owner));
  return ele;
}  // Create ()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                    seitz 12/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoNurbs27ThermoType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_nurbs27;
  Nurbs::SoNurbs27Type::setup_element_definition(definitions_nurbs27);

  std::map<std::string, Input::LineDefinition>& defs_nurbs27 = definitions_nurbs27["SONURBS27"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["NURBS27"] = defs_nurbs27["NURBS27"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          seitz 12/15 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoNurbs27ThermoType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::Nurbs::SoNurbs27, Core::FE::CellType::nurbs27>*
        actele = dynamic_cast<Discret::ELEMENTS::So3Thermo<Discret::ELEMENTS::Nurbs::SoNurbs27,
            Core::FE::CellType::nurbs27>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex20_thermo* failed");

    actele->SoNurbs27::init_jacobian_mapping(dis);
    // as an alternative we can call: So_hex27Type::Initialize(dis);
    actele->So3Thermo<Discret::ELEMENTS::Nurbs::SoNurbs27,
        Core::FE::CellType::nurbs27>::init_jacobian_mapping_special_for_nurbs(dis);
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | END nurbs27 Element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
