/*----------------------------------------------------------------------------*/
/*! \file
\brief a 2D solid-wall element with ScaTra coupling.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "wall1_scatra.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

DRT::ELEMENTS::Wall1ScatraType DRT::ELEMENTS::Wall1ScatraType::instance_;

DRT::ELEMENTS::Wall1ScatraType& DRT::ELEMENTS::Wall1ScatraType::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Wall1ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Wall1_Scatra* object = new DRT::ELEMENTS::Wall1_Scatra(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLSCATRA")
  {
    if (eledistype != "NURBS4" and eledistype != "NURBS9")
    {
      return Teuchos::rcp(new DRT::ELEMENTS::Wall1_Scatra(id, owner));
    }
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Wall1ScatraType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::Wall1_Scatra(id, owner));
}

void DRT::ELEMENTS::Wall1ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_wall;
  Wall1Type::SetupElementDefinition(definitions_wall);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["WALLSCATRA"];

  defs["QUAD4"] = defs_wall["QUAD4"];
  defs["QUAD8"] = defs_wall["QUAD8"];
  defs["QUAD9"] = defs_wall["QUAD9"];
  defs["TRI3"] = defs_wall["TRI3"];
  defs["TRI6"] = defs_wall["TRI6"];
  defs["NURBS4"] = defs_wall["NURBS4"];
  defs["NURBS9"] = defs_wall["NURBS9"];

  // add scalar transport impltype
  defs["QUAD4"].AddNamedString("TYPE");
  defs["QUAD8"].AddNamedString("TYPE");
  defs["QUAD9"].AddNamedString("TYPE");
  defs["TRI3"].AddNamedString("TYPE");
  defs["TRI6"].AddNamedString("TYPE");
  defs["NURBS4"].AddNamedString("TYPE");
  defs["NURBS9"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 01/14/|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1_Scatra::Wall1_Scatra(int id, int owner)
    : Wall1(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 01/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1_Scatra::Wall1_Scatra(const DRT::ELEMENTS::Wall1_Scatra& old)
    : Wall1(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                            vuong 01/14 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Wall1_Scatra::Clone() const
{
  DRT::ELEMENTS::Wall1_Scatra* newelement = new DRT::ELEMENTS::Wall1_Scatra(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            vuong 01/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1_Scatra::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // pack scalar transport impltype
  AddtoPack(data, impltype_);

  // add base class Element
  Wall1::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            vuong 01/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1_Scatra::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract scalar transport impltype
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(ExtractInt(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Wall1::Unpack(basedata);
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 01/14|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1_Scatra::Print(std::ostream& os) const
{
  os << "Wall1_Scatra ";
  Wall1::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Wall1_Scatra::ReadElement(
    const std::string& eletype, const std::string& eledistype, DRT::INPUT::LineDefinition* linedef)
{
  // read base element
  Wall1::ReadElement(eletype, eledistype, linedef);

  // read scalar transport implementation type
  std::string impltype;
  linedef->ExtractString("TYPE", impltype);

  if (impltype == "Undefined")
    impltype_ = INPAR::SCATRA::impltype_undefined;
  else if (impltype == "AdvReac")
    impltype_ = INPAR::SCATRA::impltype_advreac;
  else if (impltype == "CardMono")
    impltype_ = INPAR::SCATRA::impltype_cardiac_monodomain;
  else if (impltype == "Chemo")
    impltype_ = INPAR::SCATRA::impltype_chemo;
  else if (impltype == "ChemoReac")
    impltype_ = INPAR::SCATRA::impltype_chemoreac;
  else if (impltype == "Loma")
    impltype_ = INPAR::SCATRA::impltype_loma;
  else if (impltype == "RefConcReac")
    impltype_ = INPAR::SCATRA::impltype_refconcreac;
  else if (impltype == "Std")
    impltype_ = INPAR::SCATRA::impltype_std;
  else
    dserror("Invalid implementation type for Wall1_Scatra elements!");

  return true;
}
