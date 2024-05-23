/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element Type with ScaTra coupling

*----------------------------------------------------------------------*/

#include "4C_membrane_scatra_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_membrane_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatraTri3Type DRT::ELEMENTS::MembraneScatraTri3Type::instance_;

DRT::ELEMENTS::MembraneScatraTri3Type& DRT::ELEMENTS::MembraneScatraTri3Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneScatraTri3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri3>* object =
      new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatraTri3Type::Create(
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

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatraTri3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri3>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneScatraTri3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_membrane;
  MembraneTri3Type::setup_element_definition(definitions_membrane);

  std::map<std::string, INPUT::LineDefinition>& defs_membrane = definitions_membrane["MEMBRANE3"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["MEMBRANESCATRA3"];

  defs["TRI3"] =
      INPUT::LineDefinition::Builder(defs_membrane["TRI3"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  TRI 6 Element                                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatraTri6Type DRT::ELEMENTS::MembraneScatraTri6Type::instance_;

DRT::ELEMENTS::MembraneScatraTri6Type& DRT::ELEMENTS::MembraneScatraTri6Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneScatraTri6Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri6>* object =
      new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri6>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatraTri6Type::Create(
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

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatraTri6Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri6>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneScatraTri6Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_membrane;
  MembraneTri6Type::setup_element_definition(definitions_membrane);

  std::map<std::string, INPUT::LineDefinition>& defs_membrane = definitions_membrane["MEMBRANE6"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["MEMBRANESCATRA6"];

  defs["TRI6"] =
      INPUT::LineDefinition::Builder(defs_membrane["TRI6"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatraQuad4Type DRT::ELEMENTS::MembraneScatraQuad4Type::instance_;

DRT::ELEMENTS::MembraneScatraQuad4Type& DRT::ELEMENTS::MembraneScatraQuad4Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneScatraQuad4Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad4>* object =
      new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatraQuad4Type::Create(
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

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatraQuad4Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneScatraQuad4Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_membrane;
  MembraneQuad4Type::setup_element_definition(definitions_membrane);

  std::map<std::string, INPUT::LineDefinition>& defs_membrane = definitions_membrane["MEMBRANE4"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["MEMBRANESCATRA4"];

  defs["QUAD4"] =
      INPUT::LineDefinition::Builder(defs_membrane["QUAD4"]).AddNamedString("TYPE").Build();
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneScatraQuad9Type DRT::ELEMENTS::MembraneScatraQuad9Type::instance_;

DRT::ELEMENTS::MembraneScatraQuad9Type& DRT::ELEMENTS::MembraneScatraQuad9Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneScatraQuad9Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad9>* object =
      new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatraQuad9Type::Create(
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

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneScatraQuad9Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad9>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneScatraQuad9Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, INPUT::LineDefinition>> definitions_membrane;
  MembraneQuad9Type::setup_element_definition(definitions_membrane);

  std::map<std::string, INPUT::LineDefinition>& defs_membrane = definitions_membrane["MEMBRANE9"];

  std::map<std::string, INPUT::LineDefinition>& defs = definitions["MEMBRANESCATRA9"];

  defs["QUAD9"] =
      INPUT::LineDefinition::Builder(defs_membrane["QUAD9"]).AddNamedString("TYPE").Build();
}

FOUR_C_NAMESPACE_CLOSE
