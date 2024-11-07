// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_ele_immersed.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::FluidTypeImmersed Discret::Elements::FluidTypeImmersed::instance_;

Discret::Elements::FluidTypeImmersed& Discret::Elements::FluidTypeImmersed::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::FluidTypeImmersed::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::FluidImmersed* object = new Discret::Elements::FluidImmersed(-1, -1);
  object->unpack(buffer);
  return object;
}


std::shared_ptr<Core::Elements::Element> Discret::Elements::FluidTypeImmersed::create(
    const int id, const int owner)
{
  return std::make_shared<Discret::Elements::FluidImmersed>(id, owner);
}

void Discret::Elements::FluidTypeImmersed::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsimmersed = definitions["FLUIDIMMERSED"];

  defsimmersed["HEX8"] = Input::LineDefinition::Builder()
                             .add_int_vector("HEX8", 8)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            rauch 03/14|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::FluidImmersed::FluidImmersed(int id, int owner)
    : Fluid(id, owner),
      FluidImmersedBase(id, owner),
      is_immersed_(0),
      is_immersed_bdry_(0),
      has_projected_dirichletvalues_(0),
      intpoint_has_projected_divergence_(nullptr),
      stored_projected_intpoint_divergence_(nullptr)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       rauch 03/14|
 *----------------------------------------------------------------------*/
Discret::Elements::FluidImmersed::FluidImmersed(const Discret::Elements::FluidImmersed& old)
    : Fluid(old),
      FluidImmersedBase(old),
      is_immersed_(old.is_immersed_),
      is_immersed_bdry_(old.is_immersed_bdry_),
      has_projected_dirichletvalues_(old.has_projected_dirichletvalues_),
      intpoint_has_projected_divergence_(nullptr),
      stored_projected_intpoint_divergence_(nullptr)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public)  |
 |                                                          rauch 03/14 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::FluidImmersed::clone() const
{
  Discret::Elements::FluidImmersed* newelement = new Discret::Elements::FluidImmersed(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          rauch 03/14 |
 *----------------------------------------------------------------------*/
void Discret::Elements::FluidImmersed::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  Discret::Elements::Fluid::pack(data);
  // Part of immersion domain?
  add_to_pack(data, is_immersed_);
  // Part of immersion domain for immersed boundary?
  add_to_pack(data, is_immersed_bdry_);
  // has dirichletvals projected?
  add_to_pack(data, has_projected_dirichletvalues_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          rauch 03/14 |
 *----------------------------------------------------------------------*/
void Discret::Elements::FluidImmersed::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  Discret::Elements::Fluid::unpack(basedata_buffer);
  // Part of immersion domain?
  is_immersed_ = extract_int(buffer);
  // Part of immersion domain for immersed boundary?
  is_immersed_bdry_ = extract_int(buffer);
  // has dirichletvals projected?
  has_projected_dirichletvalues_ = extract_int(buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");

  return;
}

FOUR_C_NAMESPACE_CLOSE
