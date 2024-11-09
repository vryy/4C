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


Discret::Elements::SoShw6Type Discret::Elements::SoShw6Type::instance_;

Discret::Elements::SoShw6Type& Discret::Elements::SoShw6Type::instance() { return instance_; }


Core::Communication::ParObject* Discret::Elements::SoShw6Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::SoShw6(-1, -1);
  object->unpack(buffer);
  return object;
}


std::shared_ptr<Core::Elements::Element> Discret::Elements::SoShw6Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    std::shared_ptr<Core::Elements::Element> ele =
        std::make_shared<Discret::Elements::SoShw6>(id, owner);
    return ele;
  }
  return nullptr;
}


std::shared_ptr<Core::Elements::Element> Discret::Elements::SoShw6Type::create(
    const int id, const int owner)
{
  std::shared_ptr<Core::Elements::Element> ele =
      std::make_shared<Discret::Elements::SoShw6>(id, owner);
  return ele;
}


void Discret::Elements::SoShw6Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::Elements::SoShw6Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_3d_null_space(node, x0);
}

void Discret::Elements::SoShw6Type::setup_element_definition(
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
Discret::Elements::SoShw6::SoShw6(int id, int owner) : Discret::Elements::SoWeg6(id, owner)
{
  eastype_ = soshw6_easnone;
  neas_ = 0;
  optimal_parameterspace_map_ = false;
  nodes_rearranged_ = false;

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
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::SoShw6::SoShw6(const Discret::Elements::SoShw6& old)
    : Discret::Elements::SoWeg6(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::SoShw6::clone() const
{
  auto* newelement = new Discret::Elements::SoShw6(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void Discret::Elements::SoShw6::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class So_weg6 Element
  Discret::Elements::SoWeg6::pack(data);
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
void Discret::Elements::SoShw6::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class So_weg6 Element
  Discret::Elements::SoWeg6::unpack(buffer);
  // eastype_
  extract_from_pack(buffer, eastype_);
  // neas_
  extract_from_pack(buffer, neas_);
  // easdata_
  unpack_eas_data(buffer);
  // reordering
  extract_from_pack(buffer, optimal_parameterspace_map_);
  extract_from_pack(buffer, nodes_rearranged_);


  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void Discret::Elements::SoShw6::print(std::ostream& os) const
{
  os << "So_shw6 ";
  Element::print(os);
  std::cout << std::endl;
  return;
}

FOUR_C_NAMESPACE_CLOSE
