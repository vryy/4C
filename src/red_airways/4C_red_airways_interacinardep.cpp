// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_comm_pack_helpers.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_shared_ptr_from_ref.hpp"
FOUR_C_NAMESPACE_OPEN

using namespace Core::FE;

Discret::Elements::RedInterAcinarDepType Discret::Elements::RedInterAcinarDepType::instance_;

Discret::Elements::RedInterAcinarDepType& Discret::Elements::RedInterAcinarDepType::instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::Elements::RedInterAcinarDepType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::RedInterAcinarDep* object = new Discret::Elements::RedInterAcinarDep(-1, -1);
  object->unpack(buffer);
  return object;
}


/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::Elements::Element> Discret::Elements::RedInterAcinarDepType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RED_ACINAR_INTER_DEP")
  {
    std::shared_ptr<Core::Elements::Element> ele =
        std::make_shared<Discret::Elements::RedInterAcinarDep>(id, owner);
    return ele;
  }
  return nullptr;
}


/*----------------------------------------------------------------------*
 |  Create                                                              |
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::Elements::Element> Discret::Elements::RedInterAcinarDepType::create(
    const int id, const int owner)
{
  std::shared_ptr<Core::Elements::Element> ele =
      std::make_shared<Discret::Elements::RedInterAcinarDep>(id, owner);
  return ele;
}


/*----------------------------------------------------------------------*
 |  setup_element_definition                                              |
 *----------------------------------------------------------------------*/
void Discret::Elements::RedInterAcinarDepType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["RED_ACINAR_INTER_DEP"];

  defs["LINE2"] =
      Input::LineDefinition::Builder().add_int_vector("LINE2", 2).add_named_int("MAT").build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::RedInterAcinarDep::RedInterAcinarDep(int id, int owner)
    : Core::Elements::Element(id, owner)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::RedInterAcinarDep::RedInterAcinarDep(
    const Discret::Elements::RedInterAcinarDep& old)
    : Core::Elements::Element(old), elem_params_(old.elem_params_), generation_(old.generation_)
{
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedInterAcinarDep and return pointer     |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::RedInterAcinarDep::clone() const
{
  Discret::Elements::RedInterAcinarDep* newelement =
      new Discret::Elements::RedInterAcinarDep(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::Elements::RedInterAcinarDep::shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      return Core::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void Discret::Elements::RedInterAcinarDep::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Element
  Element::pack(data);

  std::map<std::string, double>::const_iterator it;

  add_to_pack(data, (elem_params_.size()));
  for (it = elem_params_.begin(); it != elem_params_.end(); it++)
  {
    add_to_pack(data, it->first);
    add_to_pack(data, it->second);
  }

  add_to_pack(data, generation_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void Discret::Elements::RedInterAcinarDep::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer base_buffer(basedata);
  Element::unpack(base_buffer);

  std::map<std::string, double> it;
  int n = 0;

  extract_from_pack(buffer, n);

  for (int i = 0; i < n; i++)
  {
    std::string name;
    double val;
    extract_from_pack(buffer, name);
    extract_from_pack(buffer, val);
    elem_params_[name] = val;
  }

  // extract generation
  extract_from_pack(buffer, generation_);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");

  return;
}



/*----------------------------------------------------------------------*
 |  Print this element (public)                             ismail 01/10|
 *----------------------------------------------------------------------*/
void Discret::Elements::RedInterAcinarDep::print(std::ostream& os) const
{
  os << "RedInterAcinarDep ";
  Element::print(os);

  return;
}


/*----------------------------------------------------------------------*
 |  Return names of visualization data                     ismail 01/10 |
 *----------------------------------------------------------------------*/
void Discret::Elements::RedInterAcinarDep::vis_names(std::map<std::string, int>& names) { return; }

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool Discret::Elements::RedInterAcinarDep::vis_data(
    const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return false;
}


/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 04/10 |
 *----------------------------------------------------------------------*/
void Discret::Elements::RedInterAcinarDep::get_params(std::string name, double& var)
{
  std::map<std::string, double>::iterator it;
  it = elem_params_.find(name);
  if (it == elem_params_.end())
  {
    FOUR_C_THROW("[%s] is not found with in the element variables", name.c_str());
    exit(1);
  }
  var = elem_params_[name];
}


/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 03/11 |
 *----------------------------------------------------------------------*/
void Discret::Elements::RedInterAcinarDep::get_params(std::string name, int& var)
{
  if (name == "Generation")
  {
    var = generation_;
  }
  else
  {
    FOUR_C_THROW("[%s] is not found with in the element INT variables", name.c_str());
    exit(1);
  }
}


/*----------------------------------------------------------------------*
 |  Get vector of lines (public)                           ismail  02/13|
 *----------------------------------------------------------------------*/
std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::RedInterAcinarDep::lines()
{
  FOUR_C_ASSERT(num_line() == 1, "RED_AIRWAY element must have one and only one line");

  return {Core::Utils::shared_ptr_from_ref(*this)};
}

FOUR_C_NAMESPACE_CLOSE
