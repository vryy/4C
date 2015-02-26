/*----------------------------------------------------------------------*/
/*!
 \file wall1_scatra.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#include "wall1_scatra.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

DRT::ELEMENTS::Wall1ScatraType DRT::ELEMENTS::Wall1ScatraType::instance_;

DRT::ELEMENTS::Wall1ScatraType& DRT::ELEMENTS::Wall1ScatraType::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::Wall1ScatraType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Wall1_Scatra* object = new DRT::ELEMENTS::Wall1_Scatra(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1ScatraType::Create( const std::string eletype,
                                                            const std::string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="WALLSCATRA" )
  {
    if ( eledistype!="NURBS4" and eledistype!="NURBS9" )
    {
      return Teuchos::rcp(new DRT::ELEMENTS::Wall1_Scatra(id,owner));
    }
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1ScatraType::Create( const int id, const int owner )
{
  return Teuchos::rcp(new DRT::ELEMENTS::Wall1_Scatra(id,owner));
}

void DRT::ELEMENTS::Wall1ScatraType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{

  std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> >  definitions_wall;
  Wall1Type::SetupElementDefinition(definitions_wall);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_wall =
      definitions_wall["WALL"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["WALLSCATRA"];

  defs["QUAD4"]=defs_wall["QUAD4"];
  defs["QUAD8"]=defs_wall["QUAD8"];
  defs["QUAD9"]=defs_wall["QUAD9"];
  defs["TRI3"]=defs_wall["TRI3"];
  defs["TRI6"]=defs_wall["TRI6"];
  defs["NURBS4"]=defs_wall["NURBS4"];
  defs["NURBS9"]=defs_wall["NURBS9"];

}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 01/08/|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1_Scatra::Wall1_Scatra(int id, int owner) :
Wall1(id,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1_Scatra::Wall1_Scatra(const DRT::ELEMENTS::Wall1_Scatra& old) :
    Wall1(old)

{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Wall1_Scatra::Clone() const
{
  DRT::ELEMENTS::Wall1_Scatra* newelement = new DRT::ELEMENTS::Wall1_Scatra(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1_Scatra::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1_Scatra::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Wall1::Unpack(basedata);
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1_Scatra::Print(std::ostream& os) const
{
  os << "Wall1_Scatra ";
  Wall1::Print(os);
  return;
}
