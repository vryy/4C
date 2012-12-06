/*!----------------------------------------------------------------------
\file So3_scatra_eletypes.cpp

<pre>
   Maintainer: Cristobal Bertoglio
               bertoglio@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "so3_scatra_eletypes.H"
#include "so3_scatra.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8ScatraType DRT::ELEMENTS::So_hex8ScatraType::instance_;


DRT::ParObject* DRT::ELEMENTS::So_hex8ScatraType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* object =
         new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ScatraType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDH8SCATRA" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>
                                                                    (id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8ScatraType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>
                                                                        (id,owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8ScatraType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_hex8;
  So_hex8Type::SetupElementDefinition(definitions_hex8);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8 =
      definitions_hex8["SOLIDH8"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["SOLIDH8SCATRA"];

  defs["HEX8"]=defs_hex8["HEX8"];
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8ScatraType::Initialize(DRT::Discretization& dis)
{

  So_hex8Type::Initialize(dis);
  return 0;
}


/*----------------------------------------------------------------------*
 |  TET 4 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet4ScatraType DRT::ELEMENTS::So_tet4ScatraType::instance_;


DRT::ParObject* DRT::ELEMENTS::So_tet4ScatraType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* object =
          new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ScatraType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDT4SCATRA" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>
                                                                    (id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4ScatraType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>
                                                                        (id,owner));
  return ele;
}

void DRT::ELEMENTS::So_tet4ScatraType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_tet4;
  So_tet4Type::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4 =
      definitions_tet4["SOLIDT4"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["SOLIDT4SCATRA"];

  defs["TET4"]=defs_tet4["TET4"];
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4ScatraType::Initialize(DRT::Discretization& dis)
{
  So_tet4Type::Initialize(dis);
  return 0;
}


/*----------------------------------------------------------------------*
 |  WEDGE 6 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_weg6ScatraType DRT::ELEMENTS::So_weg6ScatraType::instance_;


DRT::ParObject* DRT::ELEMENTS::So_weg6ScatraType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>* object =
          new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_weg6ScatraType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDW6SCATRA" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>
                                                                    (id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_weg6ScatraType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>
                                                                        (id,owner));
  return ele;
}

void DRT::ELEMENTS::So_weg6ScatraType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_weg6;
  So_weg6Type::SetupElementDefinition(definitions_weg6);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_weg6 =
      definitions_weg6["SOLIDW6"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["SOLIDW6SCATRA"];

  defs["WEDGE6"]=defs_weg6["WEDGE6"];


}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_weg6ScatraType::Initialize(DRT::Discretization& dis)
{

  So_weg6Type::Initialize(dis);
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

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStet5ScatraType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="NSTET5SCATRA" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::NStet5, DRT::Element::tet4>
                                                                    (id,owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStet5ScatraType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::NStet5, DRT::Element::tet4>
                                                                        (id,owner));
  return ele;
}

void DRT::ELEMENTS::NStet5ScatraType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
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


