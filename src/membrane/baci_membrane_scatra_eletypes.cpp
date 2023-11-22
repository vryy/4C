/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element Type with ScaTra coupling

*----------------------------------------------------------------------*/

#include "baci_membrane_scatra_eletypes.H"

#include "baci_io_linedefinition.H"
#include "baci_membrane_scatra.H"

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatra_tri3Type DRT::ELEMENTS::MembraneScatra_tri3Type::instance_;

DRT::ELEMENTS::MembraneScatra_tri3Type& DRT::ELEMENTS::MembraneScatra_tri3Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneScatra_tri3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri3>* object =
      new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_tri3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA3" && eledistype == "TRI3")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_tri3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri3>(id, owner));
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

  defs["TRI3"] =
      DRT::INPUT::LineDefinition::Builder(defs_membrane["TRI3"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  TRI 6 Element                                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatra_tri6Type DRT::ELEMENTS::MembraneScatra_tri6Type::instance_;

DRT::ELEMENTS::MembraneScatra_tri6Type& DRT::ELEMENTS::MembraneScatra_tri6Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneScatra_tri6Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri6>* object =
      new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri6>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_tri6Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA6" && eledistype == "TRI6")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri6>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_tri6Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri6>(id, owner));
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

  defs["TRI6"] =
      DRT::INPUT::LineDefinition::Builder(defs_membrane["TRI6"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatra_quad4Type DRT::ELEMENTS::MembraneScatra_quad4Type::instance_;

DRT::ELEMENTS::MembraneScatra_quad4Type& DRT::ELEMENTS::MembraneScatra_quad4Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneScatra_quad4Type::Create(
    const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad4>* object =
      new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_quad4Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA4" && eledistype == "QUAD4")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_quad4Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad4>(id, owner));
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

  defs["QUAD4"] =
      DRT::INPUT::LineDefinition::Builder(defs_membrane["QUAD4"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatra_quad9Type DRT::ELEMENTS::MembraneScatra_quad9Type::instance_;

DRT::ELEMENTS::MembraneScatra_quad9Type& DRT::ELEMENTS::MembraneScatra_quad9Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneScatra_quad9Type::Create(
    const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad9>* object =
      new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_quad9Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANESCATRA9" && eledistype == "QUAD9")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatra_quad9Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad9>(id, owner));
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

  defs["QUAD9"] =
      DRT::INPUT::LineDefinition::Builder(defs_membrane["QUAD9"]).AddNamedString("TYPE").Build();
}
