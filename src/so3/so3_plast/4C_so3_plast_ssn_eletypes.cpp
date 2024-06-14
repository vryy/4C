/*----------------------------------------------------------------------*/
/*! \file
\brief so3_plast element types
\level 2
*----------------------------------------------------------------------*/

#include "4C_so3_plast_ssn_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_plast_ssn.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *  HEX8 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of plast type                         seitz 07/13 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8PlastType Discret::ELEMENTS::SoHex8PlastType::instance_;

Discret::ELEMENTS::SoHex8PlastType& Discret::ELEMENTS::SoHex8PlastType::Instance()
{
  return instance_;
}


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoHex8PlastType::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>(id, owner));

    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8PlastType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>(id, owner));
  return ele;

}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                    seitz 07/13 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex8;
  SoHex8Type::setup_element_definition(definitions_hex8);

  std::map<std::string, Input::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = Input::LineDefinition::Builder(defs_hex8["HEX8"])
                     .add_named_string("FBAR")
                     .add_optional_named_int("NUMGP")
                     .build();

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex8PlastType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex8>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8_Plast* failed");
    // initialise all quantities
    actele->init_jacobian_mapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | END HEX8 Element
 *----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*
 *  HEX18 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
| build an instance of Plast type                         seitz 07/13 |
*----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex18PlastType Discret::ELEMENTS::SoHex18PlastType::instance_;

Discret::ELEMENTS::SoHex18PlastType& Discret::ELEMENTS::SoHex18PlastType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoHex18PlastType::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex18>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex18PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex18>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex18PlastType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex18>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex18PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX18"] = Input::LineDefinition::Builder()
                      .add_int_vector("HEX18", 18)
                      .add_named_int("MAT")
                      .add_named_string("KINEM")
                      .add_optional_named_double_vector("RAD", 3)
                      .add_optional_named_double_vector("AXI", 3)
                      .add_optional_named_double_vector("CIR", 3)
                      .add_optional_named_double_vector("FIBER1", 3)
                      .add_optional_named_double_vector("FIBER2", 3)
                      .add_optional_named_double_vector("FIBER3", 3)
                      .add_optional_named_double("STRENGTH")
                      .build();
}  // setup_element_definition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex18PlastType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex18>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex18_Plast* failed");

    actele->init_jacobian_mapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
| END HEX18 Element
*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *  HEX27 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
| build an instance of Plast type                         seitz 07/13 |
*----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex27PlastType Discret::ELEMENTS::SoHex27PlastType::instance_;

Discret::ELEMENTS::SoHex27PlastType& Discret::ELEMENTS::SoHex27PlastType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoHex27PlastType::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex27PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex27PlastType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex27>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex27PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex27;
  SoHex27Type::setup_element_definition(definitions_hex27);

  std::map<std::string, Input::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX27"] = defs_hex27["HEX27"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoHex27PlastType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex27>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex27_Plast* failed");

    actele->init_jacobian_mapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
| END HEX27 Element
*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *  tet4 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
| build an instance of Plast type                         seitz 07/13 |
*----------------------------------------------------------------------*/
Discret::ELEMENTS::SoTet4PlastType Discret::ELEMENTS::SoTet4PlastType::instance_;

Discret::ELEMENTS::SoTet4PlastType& Discret::ELEMENTS::SoTet4PlastType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoTet4PlastType::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4PlastType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoTet4PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_tet4;
  SoTet4Type::setup_element_definition(definitions_tet4);

  std::map<std::string, Input::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = defs_tet4["TET4"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoTet4PlastType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<Discret::ELEMENTS::So3Plast<Core::FE::CellType::tet4>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4_Plast* failed");

    actele->init_jacobian_mapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
| END tet4 Element
*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 *  nurbs27 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
| build an instance of Plast type                         seitz 07/13 |
*----------------------------------------------------------------------*/
Discret::ELEMENTS::SoNurbs27PlastType Discret::ELEMENTS::SoNurbs27PlastType::instance_;

Discret::ELEMENTS::SoNurbs27PlastType& Discret::ELEMENTS::SoNurbs27PlastType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::SoNurbs27PlastType::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::So3Plast<Core::FE::CellType::nurbs27>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoNurbs27PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::nurbs27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoNurbs27PlastType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::So3Plast<Core::FE::CellType::nurbs27>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoNurbs27PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["NURBS27"] = Input::LineDefinition::Builder()
                        .add_int_vector("NURBS27", 27)
                        .add_named_int("MAT")
                        .add_named_string("KINEM")
                        .build();
}  // setup_element_definition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoNurbs27PlastType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<Discret::ELEMENTS::So3Plast<Core::FE::CellType::nurbs27>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4_Plast* failed");

    actele->init_jacobian_mapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
| END nurbs27 Element
*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
