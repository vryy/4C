/*----------------------------------------------------------------------*/
/*! \file
\brief 3d TSI solid element
\level 1
\maintainer Christoph Meier

*----------------------------------------------------------------------*/

#include "so3_thermo_eletypes.H"
#include "so3_thermo.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------------*
 *  HEX8 element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 08/12 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8ThermoType DRT::ELEMENTS::So_hex8ThermoType::instance_;


/*----------------------------------------------------------------------*
 | access an instance of thermo type                                    |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8ThermoType& DRT::ELEMENTS::So_hex8ThermoType::Instance() { return instance_; }


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_hex8ThermoType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* object =
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8THERMO")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ThermoType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(id, owner));
  return ele;

}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     dano 08/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8ThermoType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex8;
  So_hex8Type::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8THERMO"];

  defs["HEX8"] = defs_hex8["HEX8"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8ThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8_thermo* failed");
    // initialise all quantities
    actele->So_hex8::InitJacobianMapping();
    // as an alternative we can call: So_hex8Type::Initialize(dis);
    actele->So3_Thermo<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>::InitJacobianMapping(dis);
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
DRT::ELEMENTS::So_hex8fbarThermoType DRT::ELEMENTS::So_hex8fbarThermoType::instance_;

DRT::ELEMENTS::So_hex8fbarThermoType& DRT::ELEMENTS::So_hex8fbarThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 05/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_hex8fbarThermoType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* object =
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 05/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8FBARTHERMO")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 05/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(id, owner));
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     dano 05/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8fbarThermoType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  // original definition
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex8fbar;

  // call setup of so3_ele
  So_hex8fbarType::SetupElementDefinition(definitions_hex8fbar);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8fbar =
      definitions_hex8fbar["SOLIDH8FBAR"];

  // templated definition
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8FBARTHERMO"];

  defs["HEX8"] = defs_hex8fbar["HEX8"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbarThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8fbar_thermo* failed");

    // initialise all quantities
    actele->So_hex8fbar::InitJacobianMapping();
    // as an alternative we can call: So_hex8fbarType::Initialize(dis);
    actele->So3_Thermo<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>::InitJacobianMapping(dis);
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
DRT::ELEMENTS::So_tet4ThermoType DRT::ELEMENTS::So_tet4ThermoType::instance_;

DRT::ELEMENTS::So_tet4ThermoType& DRT::ELEMENTS::So_tet4ThermoType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_tet4ThermoType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* object =
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT4THERMO")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 08/12 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ThermoType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(id, owner));
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | build an instance of thermo type                          dano 08/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4ThermoType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_tet4;
  So_tet4Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT4THERMO"];

  defs["TET4"] = defs_tet4["TET4"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4ThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4_thermo* failed");

    actele->So_tet4::InitJacobianMapping();
    // as an alternative we can call: So_tet4Type::Initialize(dis);
    actele->So3_Thermo<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>::InitJacobianMapping(dis);
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
DRT::ELEMENTS::So_tet10ThermoType DRT::ELEMENTS::So_tet10ThermoType::instance_;

DRT::ELEMENTS::So_tet10ThermoType& DRT::ELEMENTS::So_tet10ThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_tet10ThermoType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>* object =
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT10THERMO")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10ThermoType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>(id, owner));
  return ele;
}  // Create()


/*----------------------------------------------------------------------*
 | build an instance of thermo type                         farah 05/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10ThermoType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_tet10;
  So_tet10Type::SetupElementDefinition(definitions_tet10);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT10THERMO"];

  defs["TET10"] = defs_tet10["TET10"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          farah 05/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10ThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_tet10_thermo* failed");

    actele->So_tet10::InitJacobianMapping();
    // as an alternative we can call: So_tet4Type::Initialize(dis);
    actele->So3_Thermo<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>::InitJacobianMapping(dis);
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
DRT::ELEMENTS::So_hex27ThermoType DRT::ELEMENTS::So_hex27ThermoType::instance_;

DRT::ELEMENTS::So_hex27ThermoType& DRT::ELEMENTS::So_hex27ThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 10/13 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_hex27ThermoType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>* object =
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 10/13 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH27THERMO")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                      dano 10/13 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27ThermoType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(id, owner));
  return ele;
}  // Create ()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                     dano 10/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex27ThermoType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex27;
  So_hex27Type::SetupElementDefinition(definitions_hex27);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH27THERMO"];

  defs["HEX27"] = defs_hex27["HEX27"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                           dano 10/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex27ThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_hex27_thermo* failed");

    actele->So_hex27::InitJacobianMapping();
    // as an alternative we can call: So_hex27Type::Initialize(dis);
    actele->So3_Thermo<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>::InitJacobianMapping(dis);
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
DRT::ELEMENTS::So_hex20ThermoType DRT::ELEMENTS::So_hex20ThermoType::instance_;

DRT::ELEMENTS::So_hex20ThermoType& DRT::ELEMENTS::So_hex20ThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_hex20ThermoType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex20, DRT::Element::hex20>* object =
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex20, DRT::Element::hex20>(-1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex20ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH20THERMO")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex20, DRT::Element::hex20>(id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     farah 05/14 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex20ThermoType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex20, DRT::Element::hex20>(id, owner));
  return ele;
}  // Create ()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                    farah 05/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20ThermoType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex20;
  So_hex20Type::SetupElementDefinition(definitions_hex20);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex20 = definitions_hex20["SOLIDH20"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH20THERMO"];

  defs["HEX20"] = defs_hex20["HEX20"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          farah 05/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex20ThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex20, DRT::Element::hex20>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::So_hex20, DRT::Element::hex20>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_hex20_thermo* failed");

    actele->So_hex20::InitJacobianMapping();
    // as an alternative we can call: So_hex27Type::Initialize(dis);
    actele->So3_Thermo<DRT::ELEMENTS::So_hex20, DRT::Element::hex20>::InitJacobianMapping(dis);
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
DRT::ELEMENTS::So_nurbs27ThermoType DRT::ELEMENTS::So_nurbs27ThermoType::instance_;

DRT::ELEMENTS::So_nurbs27ThermoType& DRT::ELEMENTS::So_nurbs27ThermoType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 12/15 |
 | is called in ElementRegisterType                                     |
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::So_nurbs27ThermoType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>* object =
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>(
          -1, -1);
  object->Unpack(data);
  return object;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 12/15 |
 | is called from ParObjectFactory                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_nurbs27ThermoType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SONURBS27THERMO")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}  // Create()


/*----------------------------------------------------------------------*
 | create the new element type (public)                     seitz 12/15 |
 | virtual method of ElementType                                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_nurbs27ThermoType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>(
          id, owner));
  return ele;
}  // Create ()


/*----------------------------------------------------------------------*
 | setup the element definition (public)                    seitz 12/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_nurbs27ThermoType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_nurbs27;
  NURBS::So_nurbs27Type::SetupElementDefinition(definitions_nurbs27);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_nurbs27 =
      definitions_nurbs27["SONURBS27"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SONURBS27THERMO"];

  defs["NURBS27"] = defs_nurbs27["NURBS27"];

}  // SetupElementDefinition()


/*----------------------------------------------------------------------*
 | initialise the element (public)                          seitz 12/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_nurbs27ThermoType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;

    DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>* actele =
        dynamic_cast<
            DRT::ELEMENTS::So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_hex20_thermo* failed");

    actele->So_nurbs27::InitJacobianMapping(dis);
    // as an alternative we can call: So_hex27Type::Initialize(dis);
    actele
        ->So3_Thermo<DRT::ELEMENTS::NURBS::So_nurbs27, DRT::Element::nurbs27>::InitJacobianMapping(
            dis);
  }

  return 0;
}  // Initialize()
/*----------------------------------------------------------------------------*
 | END nurbs27 Element
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
