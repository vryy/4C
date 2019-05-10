/*!----------------------------------------------------------------------
\file membrane_scatra_eletypes.cpp

\level 3

<pre>
\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

\brief Nonlinear Membrane Finite Element Type with ScaTra coupling

*----------------------------------------------------------------------*/

#include "membrane_scatra.H"
#include "membrane_scatra_eletypes.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatra_tri3Type DRT::ELEMENTS::MembraneScatra_tri3Type::instance_;

DRT::ELEMENTS::MembraneScatra_tri3Type& DRT::ELEMENTS::MembraneScatra_tri3Type::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::MembraneScatra_tri3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<DRT::Element::tri3>* object =
      new DRT::ELEMENTS::MembraneScatra<DRT::Element::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_tri3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA3" && eledistype == "TRI3")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<DRT::Element::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_tri3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<DRT::Element::tri3>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneScatra_tri3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_membrane;
  Membrane_tri3Type::SetupElementDefinition(definitions_membrane);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_membrane =
      definitions_membrane["MEMBRANE3"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["MEMBRANESCATRA3"];

  defs["TRI3"] = defs_membrane["TRI3"];

  // add scalar transport impltype
  defs["TRI3"].AddNamedString("TYPE");

  return;
}

/*----------------------------------------------------------------------*
 |  TRI 6 Element                                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatra_tri6Type DRT::ELEMENTS::MembraneScatra_tri6Type::instance_;

DRT::ELEMENTS::MembraneScatra_tri6Type& DRT::ELEMENTS::MembraneScatra_tri6Type::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::MembraneScatra_tri6Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<DRT::Element::tri6>* object =
      new DRT::ELEMENTS::MembraneScatra<DRT::Element::tri6>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_tri6Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA6" && eledistype == "TRI6")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<DRT::Element::tri6>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_tri6Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<DRT::Element::tri6>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneScatra_tri6Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_membrane;
  Membrane_tri6Type::SetupElementDefinition(definitions_membrane);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_membrane =
      definitions_membrane["MEMBRANE6"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["MEMBRANESCATRA6"];

  defs["TRI6"] = defs_membrane["TRI6"];

  // add scalar transport impltype
  defs["TRI6"].AddNamedString("TYPE");

  return;
}

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatra_quad4Type DRT::ELEMENTS::MembraneScatra_quad4Type::instance_;

DRT::ELEMENTS::MembraneScatra_quad4Type& DRT::ELEMENTS::MembraneScatra_quad4Type::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::MembraneScatra_quad4Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<DRT::Element::quad4>* object =
      new DRT::ELEMENTS::MembraneScatra<DRT::Element::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_quad4Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA4" && eledistype == "QUAD4")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<DRT::Element::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_quad4Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<DRT::Element::quad4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneScatra_quad4Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_membrane;
  Membrane_quad4Type::SetupElementDefinition(definitions_membrane);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_membrane =
      definitions_membrane["MEMBRANE4"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["MEMBRANESCATRA4"];

  defs["QUAD4"] = defs_membrane["QUAD4"];

  // add scalar transport impltype
  defs["QUAD4"].AddNamedString("TYPE");

  return;
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatra_quad9Type DRT::ELEMENTS::MembraneScatra_quad9Type::instance_;

DRT::ELEMENTS::MembraneScatra_quad9Type& DRT::ELEMENTS::MembraneScatra_quad9Type::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::MembraneScatra_quad9Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<DRT::Element::quad9>* object =
      new DRT::ELEMENTS::MembraneScatra<DRT::Element::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_quad9Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA9" && eledistype == "QUAD9")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<DRT::Element::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_quad9Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<DRT::Element::quad9>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneScatra_quad9Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_membrane;
  Membrane_quad9Type::SetupElementDefinition(definitions_membrane);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_membrane =
      definitions_membrane["MEMBRANE9"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["MEMBRANESCATRA9"];

  defs["QUAD9"] = defs_membrane["QUAD9"];

  // add scalar transport impltype
  defs["QUAD9"].AddNamedString("TYPE");

  return;
}
