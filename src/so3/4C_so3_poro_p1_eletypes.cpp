/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the p1 (mixed) solid-poro element

 \level 2

 *----------------------------------------------------------------------*/

#include "4C_so3_poro_p1_eletypes.hpp"

#include "4C_fluid_ele_nullspace.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_poro_p1.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  HEX 8 Element                                                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8PoroP1Type Discret::ELEMENTS::SoHex8PoroP1Type::instance_;

Discret::ELEMENTS::SoHex8PoroP1Type& Discret::ELEMENTS::SoHex8PoroP1Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoHex8PoroP1Type::Create(
    const std::vector<char>& data)
{
  auto* object =
      new Discret::ELEMENTS::So3PoroP1<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3PoroP1<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8PoroP1Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3PoroP1<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoHex8PoroP1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_hex8poro;
  SoHex8PoroType::setup_element_definition(definitions_hex8poro);

  std::map<std::string, Input::LineDefinition>& defs_hex8 = definitions_hex8poro["SOLIDH8PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = defs_hex8["HEX8"];
}

void Discret::ELEMENTS::SoHex8PoroP1Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 4;
  dimns = 4;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SoHex8PoroP1Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

int Discret::ELEMENTS::SoHex8PoroP1Type::Initialize(Discret::Discretization& dis)
{
  SoHex8Type::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        Discret::ELEMENTS::So3PoroP1<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So3_Poro_P1* failed");
    actele->So3PoroP1<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>::InitElement();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  TET 4 Element                                                       |
 *----------------------------------------------------------------------*/

Discret::ELEMENTS::SoTet4PoroP1Type Discret::ELEMENTS::SoTet4PoroP1Type::instance_;

Discret::ELEMENTS::SoTet4PoroP1Type& Discret::ELEMENTS::SoTet4PoroP1Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::SoTet4PoroP1Type::Create(
    const std::vector<char>& data)
{
  auto* object =
      new Discret::ELEMENTS::So3PoroP1<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4PoroP1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
        new Discret::ELEMENTS::So3PoroP1<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
            id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoTet4PoroP1Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele = Teuchos::rcp(
      new Discret::ELEMENTS::So3PoroP1<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>(
          id, owner));
  return ele;
}

void Discret::ELEMENTS::SoTet4PoroP1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_tet4;
  SoTet4PoroType::setup_element_definition(definitions_tet4);

  std::map<std::string, Input::LineDefinition>& defs_tet4 = definitions_tet4["SOLIDT4PORO"];

  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = defs_tet4["TET4"];
}

int Discret::ELEMENTS::SoTet4PoroP1Type::Initialize(Discret::Discretization& dis)
{
  SoTet4PoroType::Initialize(dis);
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<
        Discret::ELEMENTS::So3Poro<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>*>(
        dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4_poro* failed");
    actele->So3Poro<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>::InitElement();
  }
  return 0;
}

void Discret::ELEMENTS::SoTet4PoroP1Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 4;
  dimns = 4;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SoTet4PoroP1Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

FOUR_C_NAMESPACE_CLOSE
