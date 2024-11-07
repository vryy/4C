// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_bele_vele3.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::Elements::Vele3SurfaceType Discret::Elements::Vele3SurfaceType::instance_;

Discret::Elements::Vele3SurfaceType& Discret::Elements::Vele3SurfaceType::instance()
{
  return instance_;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 05/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::Vele3Surface::Vele3Surface(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::Elements::Vele3* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}



/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
Discret::Elements::Vele3Surface::Vele3Surface(const Discret::Elements::Vele3Surface& old)
    : Core::Elements::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::Vele3Surface::clone() const
{
  Discret::Elements::Vele3Surface* newelement = new Discret::Elements::Vele3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::Elements::Vele3Surface::shape() const
{
  switch (num_node())
  {
    case 3:
      return Core::FE::CellType::tri3;
    case 4:
      return Core::FE::CellType::quad4;
    case 6:
      return Core::FE::CellType::tri6;
    case 8:
      return Core::FE::CellType::quad8;
    case 9:
      return Core::FE::CellType::quad9;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::Elements::Vele3Surface::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Vele3Surface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::Elements::Vele3Surface::unpack(Core::Communication::UnpackBuffer& buffer)
{
  FOUR_C_THROW("this Vele3Surface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::Elements::Vele3Surface::print(std::ostream& os) const
{
  os << "Vele3Surface " << Core::FE::cell_type_to_string(shape());
  Element::print(os);
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::Vele3Surface::lines()
{
  return Core::Communication::element_boundary_factory<Vele3Line, Vele3Surface>(
      Core::Communication::buildLines, *this);
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::Vele3Surface::surfaces()
{
  return {Core::Utils::shared_ptr_from_ref(*this)};
}



/*----------------------------------------------------------------------*
 |  get optimal gauss rule                                   gammi 04/07|
 *----------------------------------------------------------------------*/
Core::FE::GaussRule2D Discret::Elements::Vele3Surface::get_optimal_gaussrule(
    const Core::FE::CellType& distype) const
{
  Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::undefined;
  switch (distype)
  {
    case Core::FE::CellType::quad4:
      rule = Core::FE::GaussRule2D::quad_4point;
      break;
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
      rule = Core::FE::GaussRule2D::quad_9point;
      break;
    case Core::FE::CellType::tri3:
      rule = Core::FE::GaussRule2D::tri_3point;
      break;
    case Core::FE::CellType::tri6:
      rule = Core::FE::GaussRule2D::tri_6point;
      break;
    default:
      FOUR_C_THROW("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

FOUR_C_NAMESPACE_CLOSE
