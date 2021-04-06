/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element used for scalar transport coupling

\level 3

*/
/*---------------------------------------------------------------------------*/

#include "truss3_scatra.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3ScatraType DRT::ELEMENTS::Truss3ScatraType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3ScatraType& DRT::ELEMENTS::Truss3ScatraType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::ELEMENTS::Truss3ScatraType::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Truss3Scatra(-1, -1);
  object->Unpack(data);
  return object;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss3ScatraType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TRUSS3SCATRA")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss3Scatra(id, owner));
    return ele;
  }
  // return base class
  else
    return DRT::ELEMENTS::Truss3Type::Create(eletype, eledistype, id, owner);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss3ScatraType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss3Scatra(id, owner));
  return ele;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3ScatraType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["TRUSS3SCATRA"];

  // get definitions from standard truss element
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_truss;
  Truss3Type::SetupElementDefinition(definitions_truss);
  std::map<std::string, DRT::INPUT::LineDefinition>& defs_truss = definitions_truss["TRUSS3"];

  // copy definitions of standard truss element to truss element for scalar transport coupling
  defs["LINE2"] = defs_truss["LINE2"];

  // add scalar transport implementation type
  defs["LINE2"].AddNamedString("TYPE");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3Scatra::Truss3Scatra(int id, int owner)
    : Truss3(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Truss3Scatra::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read base element
  Truss3::ReadElement(eletype, distype, linedef);

  // read scalar transport implementation type
  std::string impltype;
  linedef->ExtractString("TYPE", impltype);

  if (impltype == "ElchDiffCond")
    impltype_ = INPAR::SCATRA::impltype_elch_diffcond;
  else
    dserror("Invalid implementation type for Truss3Scatra elements!");

  return true;
}