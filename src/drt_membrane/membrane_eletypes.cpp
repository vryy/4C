/*!----------------------------------------------------------------------
\file membrane_eletypes.cpp
\brief

\level 3

<pre>
\maintainer Fabian Bräu
            braeu@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

\brief Nonlinear Membrane Finite Element Type

*----------------------------------------------------------------------*/
#include "membrane.H"
#include "membrane_eletypes.H"

#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Membrane_tri3Type DRT::ELEMENTS::Membrane_tri3Type::instance_;

DRT::ELEMENTS::Membrane_tri3Type& DRT::ELEMENTS::Membrane_tri3Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Membrane_tri3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Membrane<DRT::Element::tri3>* object =
      new DRT::ELEMENTS::Membrane<DRT::Element::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Membrane_tri3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE3" && eledistype == "TRI3")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Membrane<DRT::Element::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Membrane_tri3Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Membrane<DRT::Element::tri3>(id, owner));
  return ele;
}

void DRT::ELEMENTS::Membrane_tri3Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

void DRT::ELEMENTS::Membrane_tri3Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Membrane_tri3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["MEMBRANE3"];

  defs["TRI3"]
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
      .AddOptionalNamedDoubleVector("FIBER3", 3);
}

/*----------------------------------------------------------------------*
 |  TRI 6 Element                                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Membrane_tri6Type DRT::ELEMENTS::Membrane_tri6Type::instance_;

DRT::ELEMENTS::Membrane_tri6Type& DRT::ELEMENTS::Membrane_tri6Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Membrane_tri6Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Membrane<DRT::Element::tri6>* object =
      new DRT::ELEMENTS::Membrane<DRT::Element::tri6>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Membrane_tri6Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE6" && eledistype == "TRI6")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Membrane<DRT::Element::tri6>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Membrane_tri6Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Membrane<DRT::Element::tri6>(id, owner));
  return ele;
}

void DRT::ELEMENTS::Membrane_tri6Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

void DRT::ELEMENTS::Membrane_tri6Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure2DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Membrane_tri6Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["MEMBRANE6"];

  defs["TRI6"]
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
      .AddOptionalNamedDoubleVector("FIBER3", 3);
}

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Membrane_quad4Type DRT::ELEMENTS::Membrane_quad4Type::instance_;

DRT::ELEMENTS::Membrane_quad4Type& DRT::ELEMENTS::Membrane_quad4Type::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::Membrane_quad4Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Membrane<DRT::Element::quad4>* object =
      new DRT::ELEMENTS::Membrane<DRT::Element::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Membrane_quad4Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE4" && eledistype == "QUAD4")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Membrane<DRT::Element::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Membrane_quad4Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Membrane<DRT::Element::quad4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::Membrane_quad4Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

void DRT::ELEMENTS::Membrane_quad4Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure2DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Membrane_quad4Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["MEMBRANE4"];

  defs["QUAD4"]
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
      .AddOptionalNamedDoubleVector("FIBER3", 3);
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Membrane_quad9Type DRT::ELEMENTS::Membrane_quad9Type::instance_;

DRT::ELEMENTS::Membrane_quad9Type& DRT::ELEMENTS::Membrane_quad9Type::Instance()
{
  return instance_;
}

DRT::ParObject* DRT::ELEMENTS::Membrane_quad9Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Membrane<DRT::Element::quad9>* object =
      new DRT::ELEMENTS::Membrane<DRT::Element::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Membrane_quad9Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE9" && eledistype == "QUAD9")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::Membrane<DRT::Element::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Membrane_quad9Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::Membrane<DRT::Element::quad9>(id, owner));
  return ele;
}

void DRT::ELEMENTS::Membrane_quad9Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

void DRT::ELEMENTS::Membrane_quad9Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure2DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Membrane_quad9Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["MEMBRANE9"];

  defs["QUAD9"]
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
      .AddOptionalNamedDoubleVector("FIBER3", 3);
}
