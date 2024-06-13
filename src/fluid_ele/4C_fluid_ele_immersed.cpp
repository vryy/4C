/*----------------------------------------------------------------------*/
/*! \file

\brief specialized immersed element used in immersed fsi

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele_immersed.hpp"

#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidTypeImmersed Discret::ELEMENTS::FluidTypeImmersed::instance_;

Discret::ELEMENTS::FluidTypeImmersed& Discret::ELEMENTS::FluidTypeImmersed::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::FluidTypeImmersed::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::FluidImmersed* object = new Discret::ELEMENTS::FluidImmersed(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::FluidTypeImmersed::Create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::FluidImmersed(id, owner));
}

void Discret::ELEMENTS::FluidTypeImmersed::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsimmersed = definitions["FLUIDIMMERSED"];

  defsimmersed["HEX8"] = Input::LineDefinition::Builder()
                             .add_int_vector("HEX8", 8)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .Build();
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            rauch 03/14|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidImmersed::FluidImmersed(int id, int owner)
    : Fluid(id, owner),
      FluidImmersedBase(id, owner),
      is_immersed_(0),
      is_immersed_bdry_(0),
      has_projected_dirichletvalues_(0),
      intpoint_has_projected_divergence_(Teuchos::null),
      stored_projected_intpoint_divergence_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       rauch 03/14|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidImmersed::FluidImmersed(const Discret::ELEMENTS::FluidImmersed& old)
    : Fluid(old),
      FluidImmersedBase(old),
      is_immersed_(old.is_immersed_),
      is_immersed_bdry_(old.is_immersed_bdry_),
      has_projected_dirichletvalues_(old.has_projected_dirichletvalues_),
      intpoint_has_projected_divergence_(Teuchos::null),
      stored_projected_intpoint_divergence_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public)  |
 |                                                          rauch 03/14 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::FluidImmersed::Clone() const
{
  Discret::ELEMENTS::FluidImmersed* newelement = new Discret::ELEMENTS::FluidImmersed(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          rauch 03/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidImmersed::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  Discret::ELEMENTS::Fluid::Pack(data);
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
void Discret::ELEMENTS::FluidImmersed::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Discret::ELEMENTS::Fluid::Unpack(basedata);
  // Part of immersion domain?
  is_immersed_ = extract_int(position, data);
  // Part of immersion domain for immersed boundary?
  is_immersed_bdry_ = extract_int(position, data);
  // has dirichletvals projected?
  has_projected_dirichletvalues_ = extract_int(position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

FOUR_C_NAMESPACE_CLOSE
