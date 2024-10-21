// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_so3_shw6.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_utils.hpp"
#include "4C_so3_weg6.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::SoShw6Type Discret::ELEMENTS::SoShw6Type::instance_;

Discret::ELEMENTS::SoShw6Type& Discret::ELEMENTS::SoShw6Type::instance() { return instance_; }


Core::Communication::ParObject* Discret::ELEMENTS::SoShw6Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::ELEMENTS::SoShw6(-1, -1);
  object->unpack(buffer);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoShw6Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::make_rcp<Discret::ELEMENTS::SoShw6>(id, owner);
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoShw6Type::create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::make_rcp<Discret::ELEMENTS::SoShw6>(id, owner);
  return ele;
}


void Discret::ELEMENTS::SoShw6Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SoShw6Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_3d_null_space(node, x0);
}

void Discret::ELEMENTS::SoShw6Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["WEDGE6"] = Input::LineDefinition::Builder()
                       .add_int_vector("WEDGE6", 6)
                       .add_named_int("MAT")
                       .add_named_string("KINEM")
                       .add_named_string("EAS")
                       .add_optional_tag("OPTORDER")
                       .add_optional_named_double_vector("RAD", 3)
                       .add_optional_named_double_vector("AXI", 3)
                       .add_optional_named_double_vector("CIR", 3)
                       .build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoShw6::SoShw6(int id, int owner) : Discret::ELEMENTS::SoWeg6(id, owner)
{
  eastype_ = soshw6_easnone;
  neas_ = 0;
  optimal_parameterspace_map_ = false;
  nodes_rearranged_ = false;

  Teuchos::RCP<const Teuchos::ParameterList> params =
      Global::Problem::instance()->get_parameter_list();
  if (params != Teuchos::null)
  {
    Discret::ELEMENTS::Utils::throw_error_fd_material_tangent(
        Global::Problem::instance()->structural_dynamic_params(), get_element_type_string());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoShw6::SoShw6(const Discret::ELEMENTS::SoShw6& old)
    : Discret::ELEMENTS::SoWeg6(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::SoShw6::clone() const
{
  auto* newelement = new Discret::ELEMENTS::SoShw6(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoShw6::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class So_weg6 Element
  Discret::ELEMENTS::SoWeg6::pack(data);
  // eastype_
  add_to_pack(data, eastype_);
  // neas_
  add_to_pack(data, neas_);
  // easdata_
  pack_eas_data(data);
  // reordering
  add_to_pack(data, optimal_parameterspace_map_);
  add_to_pack(data, nodes_rearranged_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoShw6::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class So_weg6 Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  Discret::ELEMENTS::SoWeg6::unpack(basedata_buffer);
  // eastype_
  eastype_ = static_cast<EASType>(extract_int(buffer));
  // neas_
  extract_from_pack(buffer, neas_);
  // easdata_
  unpack_eas_data(buffer);
  // reordering
  optimal_parameterspace_map_ = extract_int(buffer);
  nodes_rearranged_ = extract_int(buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoShw6::print(std::ostream& os) const
{
  os << "So_shw6 ";
  Element::print(os);
  std::cout << std::endl;
  return;
}

FOUR_C_NAMESPACE_CLOSE
