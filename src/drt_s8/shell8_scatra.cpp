/*----------------------------------------------------------------------------*/
/*! \file
\brief A 2D shell element with scatra functionality

\level 1


*/
/*---------------------------------------------------------------------------*/
#include "shell8_scatra.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::Shell8ScatraType DRT::ELEMENTS::Shell8ScatraType::instance_;


DRT::ELEMENTS::Shell8ScatraType& DRT::ELEMENTS::Shell8ScatraType::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Shell8ScatraType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Shell8_Scatra* object = new DRT::ELEMENTS::Shell8_Scatra(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 |  create pointer to shell8 scatra element (public)      schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Shell8ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SHELL8SCATRA")
  {
    return Teuchos::rcp(new DRT::ELEMENTS::Shell8_Scatra(id, owner));
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  create pointer to shell8 scatra element (public)      schmidt 09/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Shell8ScatraType::Create(const int id, const int owner)
{
  return Teuchos::rcp(new DRT::ELEMENTS::Shell8_Scatra(id, owner));
}

/*----------------------------------------------------------------------*
 |  setup of element definition (public)                  schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_shell;
  Shell8Type::SetupElementDefinition(definitions_shell);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_shell = definitions_shell["SHELL8"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SHELL8SCATRA"];

  defs["QUAD4"] = defs_shell["QUAD4"];
  defs["QUAD8"] = defs_shell["QUAD8"];
  defs["QUAD9"] = defs_shell["QUAD9"];
  defs["TRI3"] = defs_shell["TRI3"];
  defs["TRI6"] = defs_shell["TRI6"];

  // add scalar transport impltype
  defs["QUAD4"].AddNamedString("TYPE");
  defs["QUAD8"].AddNamedString("TYPE");
  defs["QUAD9"].AddNamedString("TYPE");
  defs["TRI3"].AddNamedString("TYPE");
  defs["TRI6"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8_Scatra::Shell8_Scatra(int id, int owner)
    : Shell8(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8_Scatra::Shell8_Scatra(const DRT::ELEMENTS::Shell8_Scatra& old)
    : Shell8(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Shell8_Scatra::Clone() const
{
  DRT::ELEMENTS::Shell8_Scatra* newelement = new DRT::ELEMENTS::Shell8_Scatra(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8_Scatra::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // pack impltype
  AddtoPack(data, impltype_);

  // add base class Element
  my::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data (public)                                  schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8_Scatra::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract impltype
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(ExtractInt(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  my::Unpack(basedata);
}

/*----------------------------------------------------------------------*
 |  print this element (public)                           schmidt 09/17 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8_Scatra::Print(std::ostream& os) const
{
  os << "Shell8_Scatra ";
  Element::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Shell8_Scatra::ReadElement(
    const std::string& eletype, const std::string& eledistype, DRT::INPUT::LineDefinition* linedef)
{
  // read base element
  my::ReadElement(eletype, eledistype, linedef);

  // read implementation type
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
    dserror("Invalid implementation type for Shell8_Scatra elements!");

  return true;
}
