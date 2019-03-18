/*----------------------------------------------------------------------*/
/*!
 \file so3_poro_p1_eletypes.cpp

 \brief element types of the p1 (mixed) solid-poro element

 \level 2

\maintainer  Christoph Ager
 *----------------------------------------------------------------------*/

#include "so3_poro_p1_eletypes.H"
#include "so3_poro_p1.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_utils_nullspace.H"

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                    vuong 03/12    |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8PoroP1Type DRT::ELEMENTS::So_hex8PoroP1Type::instance_;

DRT::ELEMENTS::So_hex8PoroP1Type& DRT::ELEMENTS::So_hex8PoroP1Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::So_hex8PoroP1Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* object =
      new DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8POROP1")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8PoroP1Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_hex8PoroP1Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_hex8poro;
  So_hex8PoroType::SetupElementDefinition(definitions_hex8poro);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_hex8 =
      definitions_hex8poro["SOLIDH8PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8POROP1"];

  defs["HEX8"] = defs_hex8["HEX8"];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8PoroP1Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 4;
  dimns = 4;
  nv = 3;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8PoroP1Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeFluidDNullSpace(dis, ns, x0, numdf, dimns);
}

/*----------------------------------------------------------------------*
 |  init the element (public)                        vuong 03/12        |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8PoroP1Type::Initialize(DRT::Discretization& dis)
{
  So_hex8Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So3_Poro_P1* failed");
    actele->So3_Poro_P1<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>::InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                       |
 *----------------------------------------------------------------------*/


DRT::ELEMENTS::So_tet4PoroP1Type DRT::ELEMENTS::So_tet4PoroP1Type::instance_;

DRT::ELEMENTS::So_tet4PoroP1Type& DRT::ELEMENTS::So_tet4PoroP1Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::So_tet4PoroP1Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* object =
      new DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDT4POROP1")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
        new DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4PoroP1Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(
      new DRT::ELEMENTS::So3_Poro_P1<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_tet4PoroP1Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>> definitions_tet4;
  So_tet4PoroType::SetupElementDefinition(definitions_tet4);

  std::map<std::string, DRT::INPUT::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4PORO"];

  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT4POROP1"];

  defs["TET4"] = defs_tet4["TET4"];
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4PoroP1Type::Initialize(DRT::Discretization& dis)
{
  So_tet4PoroType::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>* actele =
        dynamic_cast<DRT::ELEMENTS::So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>*>(
            dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4_poro* failed");
    actele->So3_Poro<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>::InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4PoroP1Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 4;
  dimns = 4;
  nv = 3;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4PoroP1Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeFluidDNullSpace(dis, ns, x0, numdf, dimns);
}
