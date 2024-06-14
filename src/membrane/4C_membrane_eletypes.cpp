/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 3


\brief Nonlinear Membrane Finite Element Type

*----------------------------------------------------------------------*/
#include "4C_membrane_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_membrane.hpp"
#include "4C_so3_nullspace.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneTri3Type Discret::ELEMENTS::MembraneTri3Type::instance_;

Discret::ELEMENTS::MembraneTri3Type& Discret::ELEMENTS::MembraneTri3Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::MembraneTri3Type::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Membrane<Core::FE::CellType::tri3>* object =
      new Discret::ELEMENTS::Membrane<Core::FE::CellType::tri3>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneTri3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE3" && eledistype == "TRI3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Membrane<Core::FE::CellType::tri3>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneTri3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Membrane<Core::FE::CellType::tri3>(id, owner));
  return ele;
}

void Discret::ELEMENTS::MembraneTri3Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::MembraneTri3Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void Discret::ELEMENTS::MembraneTri3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["MEMBRANE3"];

  defs["TRI3"] = Input::LineDefinition::Builder()
                     .add_int_vector("TRI3", 3)
                     .add_named_int("MAT")
                     .add_named_string("KINEM")
                     .add_named_double("THICK")
                     .add_named_string("STRESS_STRAIN")
                     .add_optional_named_double_vector("RAD", 3)
                     .add_optional_named_double_vector("AXI", 3)
                     .add_optional_named_double_vector("CIR", 3)
                     .add_optional_named_double_vector("FIBER1", 3)
                     .add_optional_named_double_vector("FIBER2", 3)
                     .add_optional_named_double_vector("FIBER3", 3)
                     .build();
}

/*----------------------------------------------------------------------*
 |  TRI 6 Element                                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneTri6Type Discret::ELEMENTS::MembraneTri6Type::instance_;

Discret::ELEMENTS::MembraneTri6Type& Discret::ELEMENTS::MembraneTri6Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::MembraneTri6Type::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Membrane<Core::FE::CellType::tri6>* object =
      new Discret::ELEMENTS::Membrane<Core::FE::CellType::tri6>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneTri6Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE6" && eledistype == "TRI6")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Membrane<Core::FE::CellType::tri6>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneTri6Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Membrane<Core::FE::CellType::tri6>(id, owner));
  return ele;
}

void Discret::ELEMENTS::MembraneTri6Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::MembraneTri6Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

void Discret::ELEMENTS::MembraneTri6Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["MEMBRANE6"];

  defs["TRI6"] = Input::LineDefinition::Builder()
                     .add_int_vector("TRI6", 6)
                     .add_named_int("MAT")
                     .add_named_string("KINEM")
                     .add_named_double("THICK")
                     .add_named_string("STRESS_STRAIN")
                     .add_optional_named_double_vector("RAD", 3)
                     .add_optional_named_double_vector("AXI", 3)
                     .add_optional_named_double_vector("CIR", 3)
                     .add_optional_named_double_vector("FIBER1", 3)
                     .add_optional_named_double_vector("FIBER2", 3)
                     .add_optional_named_double_vector("FIBER3", 3)
                     .build();
}

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneQuad4Type Discret::ELEMENTS::MembraneQuad4Type::instance_;

Discret::ELEMENTS::MembraneQuad4Type& Discret::ELEMENTS::MembraneQuad4Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::MembraneQuad4Type::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Membrane<Core::FE::CellType::quad4>* object =
      new Discret::ELEMENTS::Membrane<Core::FE::CellType::quad4>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneQuad4Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE4" && eledistype == "QUAD4")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Membrane<Core::FE::CellType::quad4>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneQuad4Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Membrane<Core::FE::CellType::quad4>(id, owner));
  return ele;
}

void Discret::ELEMENTS::MembraneQuad4Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::MembraneQuad4Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

void Discret::ELEMENTS::MembraneQuad4Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["MEMBRANE4"];

  defs["QUAD4"] = Input::LineDefinition::Builder()
                      .add_int_vector("QUAD4", 4)
                      .add_named_int("MAT")
                      .add_named_string("KINEM")
                      .add_named_double("THICK")
                      .add_named_string("STRESS_STRAIN")
                      .add_optional_named_double_vector("RAD", 3)
                      .add_optional_named_double_vector("AXI", 3)
                      .add_optional_named_double_vector("CIR", 3)
                      .add_optional_named_double_vector("FIBER1", 3)
                      .add_optional_named_double_vector("FIBER2", 3)
                      .add_optional_named_double_vector("FIBER3", 3)
                      .build();
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneQuad9Type Discret::ELEMENTS::MembraneQuad9Type::instance_;

Discret::ELEMENTS::MembraneQuad9Type& Discret::ELEMENTS::MembraneQuad9Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::MembraneQuad9Type::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Membrane<Core::FE::CellType::quad9>* object =
      new Discret::ELEMENTS::Membrane<Core::FE::CellType::quad9>(-1, -1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneQuad9Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE9" && eledistype == "QUAD9")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Membrane<Core::FE::CellType::quad9>(id, owner));
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneQuad9Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Membrane<Core::FE::CellType::quad9>(id, owner));
  return ele;
}

void Discret::ELEMENTS::MembraneQuad9Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::MembraneQuad9Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid2DNullSpace(node, x0);
}

void Discret::ELEMENTS::MembraneQuad9Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["MEMBRANE9"];

  defs["QUAD9"] = Input::LineDefinition::Builder()
                      .add_int_vector("QUAD9", 9)
                      .add_named_int("MAT")
                      .add_named_string("KINEM")
                      .add_named_double("THICK")
                      .add_named_string("STRESS_STRAIN")
                      .add_optional_named_double_vector("RAD", 3)
                      .add_optional_named_double_vector("AXI", 3)
                      .add_optional_named_double_vector("CIR", 3)
                      .add_optional_named_double_vector("FIBER1", 3)
                      .add_optional_named_double_vector("FIBER2", 3)
                      .add_optional_named_double_vector("FIBER3", 3)
                      .build();
}

FOUR_C_NAMESPACE_CLOSE
