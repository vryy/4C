/*!----------------------------------------------------------------------
\file So3_poro_eletypes.cpp

<pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "so3_poro_eletypes.H"
#include "so3_poro.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8PoroType DRT::ELEMENTS::So_hex8PoroType::instance_;


DRT::ParObject* DRT::ELEMENTS::So_hex8PoroType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* object =
        new DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDH8PORO" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>
                                                                    (id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>
                                                                        (id,owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8PoroType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_hex8;
  So_hex8Type::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8 =
      definitions_hex8["SOLIDH8"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["SOLIDH8PORO"];

  defs["HEX8"]=defs_hex8["HEX8"];
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8PoroType::Initialize(DRT::Discretization& dis)
{
  So_hex8Type::Initialize(dis);
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8> * actele =
        dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8> * >(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8_poro* failed");
    actele->So3_Poro<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>::InitElement();
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  TET 4 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet4PoroType DRT::ELEMENTS::So_tet4PoroType::instance_;


DRT::ParObject* DRT::ELEMENTS::So_tet4PoroType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* object =
        new DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDT4PORO" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>
                                                                    (id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>
                                                                   (id,owner));
  return ele;
}

void DRT::ELEMENTS::So_tet4PoroType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_tet4;
  So_tet4Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4 =
      definitions_tet4["SOLIDT4"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["SOLIDT4PORO"];

  defs["TET4"]=defs_tet4["TET4"];
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4PoroType::Initialize(DRT::Discretization& dis)
{
  So_tet4Type::Initialize(dis);
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4> * >(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4_poro* failed");
    actele->So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>::InitJacobianMapping();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  HEX 27 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_hex27PoroType DRT::ELEMENTS::So_hex27PoroType::instance_;


DRT::ParObject* DRT::ELEMENTS::So_hex27PoroType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>* object =
        new DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PoroType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDH27PORO" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>
                                                                    (id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex27PoroType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>
                                                                   (id,owner));
  return ele;
}

void DRT::ELEMENTS::So_hex27PoroType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_hex27;
  So_hex27Type::SetupElementDefinition(definitions_hex27);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex27 =
      definitions_hex27["SOLIDH27"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["SOLIDH27PORO"];

  defs["HEX27"]=defs_hex27["HEX27"];
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex27PoroType::Initialize(DRT::Discretization& dis)
{
  So_hex27PoroType::Initialize(dis);
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27> * actele =
        dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27> * >(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex27_poro* failed");
    actele->So3_Poro<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>::InitElement();
  }
  return 0;
}
