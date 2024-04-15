/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 3


\brief Nonlinear Membrane Finite Element Type

*----------------------------------------------------------------------*/
#include "baci_membrane_eletypes.hpp"

#include "baci_io_linedefinition.hpp"
#include "baci_membrane.hpp"
#include "baci_so3_nullspace.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneTri3Type DRT::ELEMENTS::MembraneTri3Type::instance_;

DRT::ELEMENTS::MembraneTri3Type& DRT::ELEMENTS::MembraneTri3Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneTri3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri3>* object =
      new DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneTri3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE3" && eledistype == "TRI3")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneTri3Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri3>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneTri3Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::MembraneTri3Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::MembraneTri3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["MEMBRANE3"];

  defs["TRI3"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("TRI3", 3)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .AddNamedDouble("THICK")
                     .AddNamedString("STRESS_STRAIN")
                     .AddOptionalNamedDoubleVector("RAD", 3)
                     .AddOptionalNamedDoubleVector("AXI", 3)
                     .AddOptionalNamedDoubleVector("CIR", 3)
                     .AddOptionalNamedDoubleVector("FIBER1", 3)
                     .AddOptionalNamedDoubleVector("FIBER2", 3)
                     .AddOptionalNamedDoubleVector("FIBER3", 3)
                     .Build();
}

/*----------------------------------------------------------------------*
 |  TRI 6 Element                                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneTri6Type DRT::ELEMENTS::MembraneTri6Type::instance_;

DRT::ELEMENTS::MembraneTri6Type& DRT::ELEMENTS::MembraneTri6Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneTri6Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri6>* object =
      new DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri6>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneTri6Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE6" && eledistype == "TRI6")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri6>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneTri6Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri6>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneTri6Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::MembraneTri6Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

void DRT::ELEMENTS::MembraneTri6Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["MEMBRANE6"];

  defs["TRI6"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("TRI6", 6)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .AddNamedDouble("THICK")
                     .AddNamedString("STRESS_STRAIN")
                     .AddOptionalNamedDoubleVector("RAD", 3)
                     .AddOptionalNamedDoubleVector("AXI", 3)
                     .AddOptionalNamedDoubleVector("CIR", 3)
                     .AddOptionalNamedDoubleVector("FIBER1", 3)
                     .AddOptionalNamedDoubleVector("FIBER2", 3)
                     .AddOptionalNamedDoubleVector("FIBER3", 3)
                     .Build();
}

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneQuad4Type DRT::ELEMENTS::MembraneQuad4Type::instance_;

DRT::ELEMENTS::MembraneQuad4Type& DRT::ELEMENTS::MembraneQuad4Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneQuad4Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad4>* object =
      new DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneQuad4Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE4" && eledistype == "QUAD4")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneQuad4Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneQuad4Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::MembraneQuad4Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

void DRT::ELEMENTS::MembraneQuad4Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["MEMBRANE4"];

  defs["QUAD4"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("QUAD4", 4)
                      .AddNamedInt("MAT")
                      .AddNamedString("KINEM")
                      .AddNamedDouble("THICK")
                      .AddNamedString("STRESS_STRAIN")
                      .AddOptionalNamedDoubleVector("RAD", 3)
                      .AddOptionalNamedDoubleVector("AXI", 3)
                      .AddOptionalNamedDoubleVector("CIR", 3)
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .AddOptionalNamedDoubleVector("FIBER2", 3)
                      .AddOptionalNamedDoubleVector("FIBER3", 3)
                      .Build();
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::MembraneQuad9Type DRT::ELEMENTS::MembraneQuad9Type::instance_;

DRT::ELEMENTS::MembraneQuad9Type& DRT::ELEMENTS::MembraneQuad9Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::MembraneQuad9Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad9>* object =
      new DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneQuad9Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE9" && eledistype == "QUAD9")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MembraneQuad9Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad9>(id, owner));
  return ele;
}

void DRT::ELEMENTS::MembraneQuad9Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::MembraneQuad9Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

void DRT::ELEMENTS::MembraneQuad9Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["MEMBRANE9"];

  defs["QUAD9"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("QUAD9", 9)
                      .AddNamedInt("MAT")
                      .AddNamedString("KINEM")
                      .AddNamedDouble("THICK")
                      .AddNamedString("STRESS_STRAIN")
                      .AddOptionalNamedDoubleVector("RAD", 3)
                      .AddOptionalNamedDoubleVector("AXI", 3)
                      .AddOptionalNamedDoubleVector("CIR", 3)
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .AddOptionalNamedDoubleVector("FIBER2", 3)
                      .AddOptionalNamedDoubleVector("FIBER3", 3)
                      .Build();
}

FOUR_C_NAMESPACE_CLOSE
