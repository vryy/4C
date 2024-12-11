// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_so3_sh18.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_so3_utils.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::SoSh18Type Discret::Elements::SoSh18Type::instance_;

Discret::Elements::SoSh18Type& Discret::Elements::SoSh18Type::instance() { return instance_; }
namespace
{
  const std::string name = Discret::Elements::SoSh18Type::instance().name();
}

Core::Communication::ParObject* Discret::Elements::SoSh18Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::SoSh18(-1, -1);
  object->unpack(buffer);
  return object;
}


std::shared_ptr<Core::Elements::Element> Discret::Elements::SoSh18Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    std::shared_ptr<Core::Elements::Element> ele =
        std::make_shared<Discret::Elements::SoSh18>(id, owner);
    return ele;
  }

  return nullptr;
}


std::shared_ptr<Core::Elements::Element> Discret::Elements::SoSh18Type::create(
    const int id, const int owner)
{
  std::shared_ptr<Core::Elements::Element> ele =
      std::make_shared<Discret::Elements::SoSh18>(id, owner);
  return ele;
}

void Discret::Elements::SoSh18Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX18"] = Input::LineDefinition::Builder()
                      .add_named_int_vector("HEX18", 18)
                      .add_named_int("MAT")
                      .add_named_string("KINEM")
                      .add_named_string("TSL")
                      .add_named_string("MEL")
                      .add_named_string("CTL")
                      .add_named_string("VOL")
                      .add_optional_named_double_vector("RAD", 3)
                      .add_optional_named_double_vector("AXI", 3)
                      .add_optional_named_double_vector("CIR", 3)
                      .add_optional_named_double_vector("FIBER1", 3)
                      .add_optional_named_double_vector("FIBER2", 3)
                      .add_optional_named_double_vector("FIBER3", 3)
                      .add_optional_named_double("STRENGTH")
                      .build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           seitz 11/14 |
 *----------------------------------------------------------------------*/
Discret::Elements::SoSh18::SoSh18(int id, int owner) : SoBase(id, owner), SoHex18(id, owner)
{
  std::shared_ptr<const Teuchos::ParameterList> params =
      Global::Problem::instance()->get_parameter_list();
  if (params != nullptr)
  {
    Discret::Elements::Utils::throw_error_fd_material_tangent(
        Global::Problem::instance()->structural_dynamic_params(), get_element_type_string());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      seitz 11/14 |
 *----------------------------------------------------------------------*/
Discret::Elements::SoSh18::SoSh18(const Discret::Elements::SoSh18& old)
    : SoBase(old),
      SoHex18(old),
      dsg_shear_(old.dsg_shear_),
      dsg_membrane_(old.dsg_membrane_),
      dsg_ctl_(old.dsg_ctl_),
      eas_(old.eas_)
{
  setup_dsg();
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::SoSh18::clone() const
{
  auto* newelement = new Discret::Elements::SoSh18(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::Elements::SoSh18::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  SoBase::pack(data);

  // detJ_
  add_to_pack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i) add_to_pack(data, invJ_[i]);

  // element technology bools
  add_to_pack(data, dsg_shear_);
  add_to_pack(data, dsg_membrane_);
  add_to_pack(data, dsg_ctl_);
  add_to_pack(data, eas_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::Elements::SoSh18::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  SoBase::unpack(buffer);

  // detJ_
  extract_from_pack(buffer, detJ_);
  // invJ_
  int size = 0;
  extract_from_pack(buffer, size);
  invJ_.resize(size, Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>(true));
  for (int i = 0; i < size; ++i) extract_from_pack(buffer, invJ_[i]);

  // element technology bools
  extract_from_pack(buffer, dsg_shear_);
  extract_from_pack(buffer, dsg_membrane_);
  extract_from_pack(buffer, dsg_ctl_);
  extract_from_pack(buffer, eas_);
  setup_dsg();


  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::Elements::SoSh18::print(std::ostream& os) const
{
  os << "So_sh18 ";
  Element::print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
