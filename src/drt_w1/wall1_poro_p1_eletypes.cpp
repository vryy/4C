/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro_p1_eletypes.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#include "wall1_poro_p1_eletypes.H"
#include "wall1_poro_p1.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                          vuong 07/13 |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallQuad4PoroP1Type DRT::ELEMENTS::WallQuad4PoroP1Type::instance_;


DRT::ParObject* DRT::ELEMENTS::WallQuad4PoroP1Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>* object = new DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>(-1,-1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 *                                                           vuong 07/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad4PoroP1Type::Create( const std::string eletype,
                                                            const std::string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="WALLQ4POROP1" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>(id,owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *                                                           vuong 07/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad4PoroP1Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>(id,owner));
  return ele;
}

/*----------------------------------------------------------------------*
 *                                                           vuong 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallQuad4PoroP1Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_wallporo;
  WallQuad4PoroType::SetupElementDefinition(definitions_wallporo);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_wallporo =
      definitions_wallporo["WALLQ4PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["WALLQ4POROP1"];

  defs["QUAD4"]=defs_wallporo["QUAD4"];
}

/*----------------------------------------------------------------------*
 |  init the element (public)                              vuong 07/13     |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::WallQuad4PoroP1Type::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>* actele = dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>*>(dis.lColElement(i));
    if (!actele) dserror("cast to Wall1_PoroP1* failed");
    actele->InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                       vuong 07/13 |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::WallQuad9PoroP1Type DRT::ELEMENTS::WallQuad9PoroP1Type::instance_;


DRT::ParObject* DRT::ELEMENTS::WallQuad9PoroP1Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>* object = new DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>(-1,-1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 *                                                           vuong 07/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad9PoroP1Type::Create( const std::string eletype,
                                                            const std::string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="WALLQ9POROP1" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>(id,owner));
    return ele;
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *                                                           vuong 07/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::WallQuad9PoroP1Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>(id,owner));
  return ele;
}

/*----------------------------------------------------------------------*
 *                                                           vuong 07/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::WallQuad9PoroP1Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_wallporo;
  WallQuad9PoroType::SetupElementDefinition(definitions_wallporo);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_wallporo =
      definitions_wallporo["WALLQ9PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["WALLQ9POROP1"];

  defs["QUAD9"]=defs_wallporo["QUAD9"];
}

/*----------------------------------------------------------------------*
 |  init the element (public)                             vuong 07/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::WallQuad9PoroP1Type::Initialize(DRT::Discretization& dis)
{
  DRT::ELEMENTS::Wall1Type::Initialize(dis);
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>* actele = dynamic_cast<DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>*>(dis.lColElement(i));
    if (!actele) dserror("cast to Wall1_PoroP1* failed");
    actele->InitElement();
  }
  return 0;
}
