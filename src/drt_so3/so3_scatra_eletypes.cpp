/*!----------------------------------------------------------------------
\file so3_scatra_eletypes.cpp

\brief Solid-scatra elements types

\level 2

<pre>
   \maintainer Thon Moritz
               thon@mhpc.mw.tum.de
               089 - 289-10264
</pre>

*----------------------------------------------------------------------*/

#include "so3_scatra.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8ScatraType DRT::ELEMENTS::So_hex8ScatraType::instance_;

DRT::ELEMENTS::So_hex8ScatraType& DRT::ELEMENTS::So_hex8ScatraType::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::So_hex8ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex8;
  So_hex8Type::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8SCATRA"];

  defs["HEX8"] = defs_hex8["HEX8"];

  // add scalar transport ImplType
  defs["HEX8"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8ScatraType::Initialize(DRT::Discretization& dis)
{
  So_hex8Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8_scatra* failed");
    actele->InitElement();
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  HEX 8 fbar Element                                       Thon 12/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8fbarScatraType DRT::ELEMENTS::So_hex8fbarScatraType::instance_;

DRT::ELEMENTS::So_hex8fbarScatraType& DRT::ELEMENTS::So_hex8fbarScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_hex8fbarScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8FBARSCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8fbarScatraType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8fbarScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex8;
  So_hex8fbarType::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8 = definitions_hex8["SOLIDH8FBAR"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8FBARSCATRA"];

  defs["HEX8"] = defs_hex8["HEX8"];

  // add scalar transport ImplType
  defs["HEX8"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbarScatraType::Initialize(DRT::Discretization& dis)
{
  So_hex8fbarType::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8fbar_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  HEX 27 Solid Scatra Element                              thon 12/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex27ScatraType DRT::ELEMENTS::So_hex27ScatraType::instance_;

DRT::ELEMENTS::So_hex27ScatraType& DRT::ELEMENTS::So_hex27ScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_hex27ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH27SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex27ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex27;
  So_hex27Type::SetupElementDefinition(definitions_hex27);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex27 = definitions_hex27["SOLIDH27"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH27SCATRA"];

  defs["HEX27"] = defs_hex27["HEX27"];

  // add scalar transport ImplType
  defs["HEX27"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex27ScatraType::Initialize(DRT::Discretization& dis)
{
  So_hex27Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_hex27_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet4ScatraType DRT::ELEMENTS::So_tet4ScatraType::instance_;

DRT::ELEMENTS::So_tet4ScatraType& DRT::ELEMENTS::So_tet4ScatraType::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::So_tet4ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT4SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet4ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_tet4;
  So_tet4Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT4SCATRA"];

  defs["TET4"] = defs_tet4["TET4"];

  // add scalar transport ImplType
  defs["TET4"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4ScatraType::Initialize(DRT::Discretization& dis)
{
  So_tet4Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 10 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet10ScatraType DRT::ELEMENTS::So_tet10ScatraType::instance_;

DRT::ELEMENTS::So_tet10ScatraType& DRT::ELEMENTS::So_tet10ScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::So_tet10ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT10SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet10ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_tet10;
  So_tet10Type::SetupElementDefinition(definitions_tet10);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet10 = definitions_tet10["SOLIDT10"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT10SCATRA"];

  defs["TET10"] = defs_tet10["TET10"];

  // add scalar transport ImplType
  defs["TET10"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10ScatraType::Initialize(DRT::Discretization& dis)
{
  So_tet10Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_tet10_scatra* failed");
    actele->InitElement();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  WEDGE 6 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_weg6ScatraType DRT::ELEMENTS::So_weg6ScatraType::instance_;

DRT::ELEMENTS::So_weg6ScatraType& DRT::ELEMENTS::So_weg6ScatraType::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::So_weg6ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>* object =
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_weg6ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDW6SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_weg6ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_weg6ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_weg6;
  So_weg6Type::SetupElementDefinition(definitions_weg6);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_weg6 = definitions_weg6["SOLIDW6"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDW6SCATRA"];

  defs["WEDGE6"] = defs_weg6["WEDGE6"];

  // add scalar transport ImplType
  defs["WEDGE6"].AddNamedString("TYPE");
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_weg6ScatraType::Initialize(DRT::Discretization& dis)
{
  So_weg6Type::Initialize(dis);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_weg6_scatra* failed");
    actele->InitElement();
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  NSTET 5 Element                                       |
 *----------------------------------------------------------------------*/

/*
DRT::ELEMENTS::NStet5ScatraType DRT::ELEMENTS::NStet5ScatraType::instance_;


DRT::ParObject* DRT::ELEMENTS::NStet5ScatraType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::NStet5, DRT::Element::tet4>* object =
          new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::NStet5, DRT::Element::tet4>(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStet5ScatraType::Create( const std::string eletype,
                                                            const std::string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="NSTET5SCATRA" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new
DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::NStet5, DRT::Element::tet4> (id,owner)); return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStet5ScatraType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::NStet5,
DRT::Element::tet4> (id,owner)); return ele;
}

void DRT::ELEMENTS::NStet5ScatraType::SetupElementDefinition(
std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_nstet5;
  NStet5Type::SetupElementDefinition(definitions_nstet5);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_nstet5 =
      definitions_nstet5["NSTET5"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["NSTET5SCATRA"];

  defs["TET4"]=defs_nstet5["TET4"];


}

*/

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
/*

int DRT::ELEMENTS::NStet5ScatraType::Initialize(DRT::Discretization& dis)
{
  NStet5Type::Initialize(dis);
  return 0;
}

*/
