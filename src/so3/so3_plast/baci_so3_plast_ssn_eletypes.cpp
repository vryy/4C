/*----------------------------------------------------------------------*/
/*! \file
\brief so3_plast element types
\level 2
*----------------------------------------------------------------------*/

#include "baci_so3_plast_ssn_eletypes.H"

#include "baci_io_linedefinition.H"
#include "baci_so3_plast_ssn.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *  HEX8 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of plast type                         seitz 07/13 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8PlastType DRT::ELEMENTS::So_hex8PlastType::instance_;

DRT::ELEMENTS::So_hex8PlastType& DRT::ELEMENTS::So_hex8PlastType::Instance() { return instance_; }


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::So_hex8PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex8>(id, owner));

    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 07/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex8>(id, owner));
  return ele;

}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                    seitz 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8PlastType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex8;
  So_hex8Type::SetupElementDefinition(definitions_hex8);

  std::map<std::string, INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX8"] = INPUT::LineDefinition::Builder(defs_hex8["HEX8"])
                     .AddNamedString("FBAR")
                     .AddOptionalNamedInt("NUMGP")
                     .Build();

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex8>*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8_Plast* failed");
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
DRT::ELEMENTS::So_hex18PlastType DRT::ELEMENTS::So_hex18PlastType::instance_;

DRT::ELEMENTS::So_hex18PlastType& DRT::ELEMENTS::So_hex18PlastType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::So_hex18PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex18>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex18PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex18>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex18PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex18>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex18PlastType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX18"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("HEX18", 18)
                      .AddNamedInt("MAT")
                      .AddNamedString("KINEM")
                      .AddOptionalNamedDoubleVector("RAD", 3)
                      .AddOptionalNamedDoubleVector("AXI", 3)
                      .AddOptionalNamedDoubleVector("CIR", 3)
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .AddOptionalNamedDoubleVector("FIBER2", 3)
                      .AddOptionalNamedDoubleVector("FIBER3", 3)
                      .AddOptionalNamedDouble("STRENGTH")
                      .Build();
}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex18PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex18>*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex18_Plast* failed");

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
DRT::ELEMENTS::So_hex27PlastType DRT::ELEMENTS::So_hex27PlastType::instance_;

DRT::ELEMENTS::So_hex27PlastType& DRT::ELEMENTS::So_hex27PlastType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::So_hex27PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex27>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex27PlastType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_hex27;
  So_hex27Type::SetupElementDefinition(definitions_hex27);

  std::map<std::string, INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX27"] = defs_hex27["HEX27"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex27PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::hex27>*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex27_Plast* failed");

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
DRT::ELEMENTS::So_tet4PlastType DRT::ELEMENTS::So_tet4PlastType::instance_;

DRT::ELEMENTS::So_tet4PlastType& DRT::ELEMENTS::So_tet4PlastType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::So_tet4PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::tet4>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4PlastType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_tet4;
  So_tet4Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["TET4"] = defs_tet4["TET4"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::tet4>*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4_Plast* failed");

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
DRT::ELEMENTS::So_nurbs27PlastType DRT::ELEMENTS::So_nurbs27PlastType::instance_;

DRT::ELEMENTS::So_nurbs27PlastType& DRT::ELEMENTS::So_nurbs27PlastType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called in ElementRegisterType                                     |
*----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::So_nurbs27PlastType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::nurbs27>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| is called from ParObjectFactory                                      |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_nurbs27PlastType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::nurbs27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
| create the new element type (public)                     seitz 07/13 |
| virtual method of ElementType                                        |
*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_nurbs27PlastType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::nurbs27>(id, owner));
  return ele;
}  // Create()


/*---------------------------------------------------------------------*
|                                                          seitz 07/13 |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_nurbs27PlastType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["NURBS27"] = INPUT::LineDefinition::Builder()
                        .AddIntVector("NURBS27", 27)
                        .AddNamedInt("MAT")
                        .AddNamedString("KINEM")
                        .Build();
}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
| initialise the element (public)                          seitz 07/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_nurbs27PlastType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    auto* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Plast<CORE::FE::CellType::nurbs27>*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4_Plast* failed");

    actele->InitJacobianMapping();
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
| END nurbs27 Element
*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/

BACI_NAMESPACE_CLOSE
