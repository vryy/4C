// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_membrane_eletypes.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_membrane.hpp"
#include "4C_so3_nullspace.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
Discret::Elements::MembraneTri3Type Discret::Elements::MembraneTri3Type::instance_;

Discret::Elements::MembraneTri3Type& Discret::Elements::MembraneTri3Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::MembraneTri3Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::Membrane<Core::FE::CellType::tri3>* object =
      new Discret::Elements::Membrane<Core::FE::CellType::tri3>(-1, -1);
  object->unpack(buffer);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneTri3Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE3" && eledistype == "TRI3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::Elements::Membrane<Core::FE::CellType::tri3>>(id, owner);
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneTri3Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::Elements::Membrane<Core::FE::CellType::tri3>>(id, owner);
  return ele;
}

void Discret::Elements::MembraneTri3Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

Core::LinAlg::SerialDenseMatrix Discret::Elements::MembraneTri3Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_3d_null_space(node, x0);
}

void Discret::Elements::MembraneTri3Type::setup_element_definition(
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
Discret::Elements::MembraneTri6Type Discret::Elements::MembraneTri6Type::instance_;

Discret::Elements::MembraneTri6Type& Discret::Elements::MembraneTri6Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::MembraneTri6Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::Membrane<Core::FE::CellType::tri6>* object =
      new Discret::Elements::Membrane<Core::FE::CellType::tri6>(-1, -1);
  object->unpack(buffer);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneTri6Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE6" && eledistype == "TRI6")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::Elements::Membrane<Core::FE::CellType::tri6>>(id, owner);
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneTri6Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::Elements::Membrane<Core::FE::CellType::tri6>>(id, owner);
  return ele;
}

void Discret::Elements::MembraneTri6Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

Core::LinAlg::SerialDenseMatrix Discret::Elements::MembraneTri6Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_2d_null_space(node, x0);
}

void Discret::Elements::MembraneTri6Type::setup_element_definition(
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
Discret::Elements::MembraneQuad4Type Discret::Elements::MembraneQuad4Type::instance_;

Discret::Elements::MembraneQuad4Type& Discret::Elements::MembraneQuad4Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::MembraneQuad4Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::Membrane<Core::FE::CellType::quad4>* object =
      new Discret::Elements::Membrane<Core::FE::CellType::quad4>(-1, -1);
  object->unpack(buffer);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneQuad4Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE4" && eledistype == "QUAD4")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::Elements::Membrane<Core::FE::CellType::quad4>>(id, owner);
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneQuad4Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::Elements::Membrane<Core::FE::CellType::quad4>>(id, owner);
  return ele;
}

void Discret::Elements::MembraneQuad4Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

Core::LinAlg::SerialDenseMatrix Discret::Elements::MembraneQuad4Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_2d_null_space(node, x0);
}

void Discret::Elements::MembraneQuad4Type::setup_element_definition(
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
Discret::Elements::MembraneQuad9Type Discret::Elements::MembraneQuad9Type::instance_;

Discret::Elements::MembraneQuad9Type& Discret::Elements::MembraneQuad9Type::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::MembraneQuad9Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::Membrane<Core::FE::CellType::quad9>* object =
      new Discret::Elements::Membrane<Core::FE::CellType::quad9>(-1, -1);
  object->unpack(buffer);
  return object;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneQuad9Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "MEMBRANE9" && eledistype == "QUAD9")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::Elements::Membrane<Core::FE::CellType::quad9>>(id, owner);
    return ele;
  }
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::Elements::MembraneQuad9Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::Elements::Membrane<Core::FE::CellType::quad9>>(id, owner);
  return ele;
}

void Discret::Elements::MembraneQuad9Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;  // number of degrees of freedom per node
  dimns = 3;  // number of nullspace vectors
  nv = 3;     // default value for no. of velocity dofs
}

Core::LinAlg::SerialDenseMatrix Discret::Elements::MembraneQuad9Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_2d_null_space(node, x0);
}

void Discret::Elements::MembraneQuad9Type::setup_element_definition(
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
