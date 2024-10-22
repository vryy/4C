// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_membrane.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  LINE 2 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneLine2Type Discret::ELEMENTS::MembraneLine2Type::instance_;

Discret::ELEMENTS::MembraneLine2Type& Discret::ELEMENTS::MembraneLine2Type::instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |  LINE 3 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::MembraneLine3Type Discret::ELEMENTS::MembraneLine3Type::instance_;

Discret::ELEMENTS::MembraneLine3Type& Discret::ELEMENTS::MembraneLine3Type::instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 |  constructor (public)                                   fbraeu 06/16 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::MembraneLine<distype>::MembraneLine(int id, int owner, int nnode,
    const int* nodeids, Core::Nodes::Node** nodes, Discret::ELEMENTS::Membrane<distype>* parent,
    const int lline)
    : Core::Elements::FaceElement(id, owner), intpointsline_(Core::FE::GaussRule1D::line_2point)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lline);
  switch (Core::FE::DisTypeToFaceShapeType<distype>::shape)
  {
    case Core::FE::CellType::line2:
    {
      Core::FE::GaussRule1D gaussrule = Core::FE::GaussRule1D::line_2point;
      // get gauss integration points
      intpointsline_ = Core::FE::IntegrationPoints1D(gaussrule);
      break;
    }
    case Core::FE::CellType::line3:
    {
      Core::FE::GaussRule1D gaussrule = Core::FE::GaussRule1D::line_3point;
      // get gauss integration points
      intpointsline_ = Core::FE::IntegrationPoints1D(gaussrule);
      break;
    }
    break;
    default:
      FOUR_C_THROW("shape type unknown!\n");
      break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  copy-constructor (public)                              fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::MembraneLine<distype>::MembraneLine(
    const Discret::ELEMENTS::MembraneLine<distype>& old)
    : Core::Elements::FaceElement(old), intpointsline_(old.intpointsline_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::MembraneLine<distype>::clone() const
{
  Discret::ELEMENTS::MembraneLine<distype>* newelement =
      new Discret::ELEMENTS::MembraneLine<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::FE::CellType Discret::ELEMENTS::MembraneLine<distype>::shape() const
{
  return Core::FE::DisTypeToFaceShapeType<distype>::shape;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             fb 09/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::MembraneLine<distype>::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this membrane line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                        fbraeu 06/165 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::MembraneLine<distype>::unpack(Core::Communication::UnpackBuffer& buffer)
{
  FOUR_C_THROW("this membrane line element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             fbraeu 06/16|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::MembraneLine<distype>::print(std::ostream& os) const
{
  os << "MembraneLine ";
  os << " discretization type: "
     << Core::FE::cell_type_to_string(Core::FE::DisTypeToFaceShapeType<distype>::shape).c_str();
  Element::print(os);
  return;
}

template class Discret::ELEMENTS::MembraneLine<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::MembraneLine<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::MembraneLine<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::MembraneLine<Core::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
