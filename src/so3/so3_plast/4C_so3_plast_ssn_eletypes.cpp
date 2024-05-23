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
DRT::ELEMENTS::SoHex8PlastType DRT::ELEMENTS::SoHex8PlastType::instance_;

DRT::ELEMENTS::SoHex8PlastType& DRT::ELEMENTS::SoHex8PlastType::Instance() { return instance_; }


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::SoHex8PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex8PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex8>(id, owner));

    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex8PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex8>(id, owner));
  return ele;

}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                    seitz 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8;
  SoHex8Type::setup_element_definition(definitions_hex8);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"])
                     .AddNamedString("FBAR")
                     .AddOptionalNamedInt("NUMGP")
                     .Build();

}  // setup_element_definition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex8PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex8>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex8_Plast* failed");
    // initialise all quantities
    actele->InitJacobianMapping();
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
DRT::ELEMENTS::SoHex18PlastType DRT::ELEMENTS::SoHex18PlastType::instance_;

DRT::ELEMENTS::SoHex18PlastType& DRT::ELEMENTS::SoHex18PlastType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::SoHex18PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex18>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex18PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex18>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex18PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex18>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex18PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX18"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("HEX18", 18)
                      .AddNamedInt("MAT")
                      .AddNamedString("KINEM")
                      .add_optional_named_double_vector("RAD", 3)
                      .add_optional_named_double_vector("AXI", 3)
                      .add_optional_named_double_vector("CIR", 3)
                      .add_optional_named_double_vector("FIBER1", 3)
                      .add_optional_named_double_vector("FIBER2", 3)
                      .add_optional_named_double_vector("FIBER3", 3)
                      .add_optional_named_double("STRENGTH")
                      .Build();
}  // setup_element_definition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex18PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex18>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex18_Plast* failed");

    actele->InitJacobianMapping();
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
DRT::ELEMENTS::SoHex27PlastType DRT::ELEMENTS::SoHex27PlastType::instance_;

DRT::ELEMENTS::SoHex27PlastType& DRT::ELEMENTS::SoHex27PlastType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::SoHex27PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex27PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoHex27PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex27>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex27PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex27;
  SoHex27Type::setup_element_definition(definitions_hex27);

  std::map<std::string, INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX27"] = defs_hex27["HEX27"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoHex27PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Plast<CORE::FE::CellType::hex27>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex27_Plast* failed");

    actele->InitJacobianMapping();
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
DRT::ELEMENTS::SoTet4PlastType DRT::ELEMENTS::SoTet4PlastType::instance_;

DRT::ELEMENTS::SoTet4PlastType& DRT::ELEMENTS::SoTet4PlastType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::SoTet4PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoTet4PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoTet4PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::tet4>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoTet4PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet4;
  SoTet4Type::setup_element_definition(definitions_tet4);

  std::map<std::string, INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = defs_tet4["TET4"];

}  // setup_element_definition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoTet4PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Plast<CORE::FE::CellType::tet4>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4_Plast* failed");

    actele->InitJacobianMapping();
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
DRT::ELEMENTS::SoNurbs27PlastType DRT::ELEMENTS::SoNurbs27PlastType::instance_;

DRT::ELEMENTS::SoNurbs27PlastType& DRT::ELEMENTS::SoNurbs27PlastType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::SoNurbs27PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::nurbs27>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoNurbs27PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::nurbs27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoNurbs27PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3Plast<CORE::FE::CellType::nurbs27>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoNurbs27PlastType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["NURBS27"] = INPUT::LineDefinition::Builder()
                        .AddIntVector("NURBS27", 27)
                        .AddNamedInt("MAT")
                        .AddNamedString("KINEM")
                        .Build();
}  // setup_element_definition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoNurbs27PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3Plast<CORE::FE::CellType::nurbs27>*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4_Plast* failed");

    actele->InitJacobianMapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
| END nurbs27 Element
*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
