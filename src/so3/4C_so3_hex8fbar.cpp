// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_so3_hex8fbar.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_prestress.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::SoHex8fbarType Discret::Elements::SoHex8fbarType::instance_;

Discret::Elements::SoHex8fbarType& Discret::Elements::SoHex8fbarType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::SoHex8fbarType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::SoHex8fbar(-1, -1);
  object->unpack(buffer);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::Elements::SoHex8fbarType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::Elements::SoHex8fbar>(id, owner);
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::Elements::SoHex8fbarType::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::Elements::SoHex8fbar>(id, owner);
  return ele;
}


void Discret::Elements::SoHex8fbarType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
  np = 0;
}

Core::LinAlg::SerialDenseMatrix Discret::Elements::SoHex8fbarType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_3d_null_space(node, x0);
}

void Discret::Elements::SoHex8fbarType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = Input::LineDefinition::Builder()
                     .add_int_vector("HEX8", 8)
                     .add_named_int("MAT")
                     .add_named_string("KINEM")
                     .add_optional_named_double_vector("RAD", 3)
                     .add_optional_named_double_vector("AXI", 3)
                     .add_optional_named_double_vector("CIR", 3)
                     .add_optional_named_double_vector("FIBER1", 3)
                     .add_optional_named_double_vector("FIBER2", 3)
                     .add_optional_named_double_vector("FIBER3", 3)
                     .add_optional_named_double("GROWTHTRIG")
                     .build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 07/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::SoHex8fbar::SoHex8fbar(int id, int owner) : Discret::Elements::SoHex8(id, owner)
{
  if (Prestress::is_mulf(pstype_))
    prestress_ = Teuchos::make_rcp<Discret::Elements::PreStress>(NUMNOD_SOH8, NUMGPT_SOH8 + 1);

  Teuchos::RCP<const Teuchos::ParameterList> params =
      Global::Problem::instance()->get_parameter_list();
  if (params != Teuchos::null)
  {
    Discret::Elements::Utils::throw_error_fd_material_tangent(
        Global::Problem::instance()->structural_dynamic_params(), get_element_type_string());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        popp 07/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::SoHex8fbar::SoHex8fbar(const Discret::Elements::SoHex8fbar& old)
    : Discret::Elements::SoHex8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            popp 07/10|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::SoHex8fbar::clone() const
{
  auto* newelement = new Discret::Elements::SoHex8fbar(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            popp 07/10|
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex8fbar::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class So_hex8 Element
  Discret::Elements::SoHex8::pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            popp 07/10|
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex8fbar::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class So_hex8 Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  Discret::Elements::SoHex8::unpack(basedata_buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                               popp 07/10|
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex8fbar::print(std::ostream& os) const
{
  os << "So_hex8fbar ";
  Element::print(os);
  std::cout << std::endl;
  return;
}

FOUR_C_NAMESPACE_CLOSE
