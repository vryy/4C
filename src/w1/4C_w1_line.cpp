/*----------------------------------------------------------------------------*/
/*! \file
\brief Line definitions for wall1 element.

\level 1


*/
/*---------------------------------------------------------------------------*/

#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_w1.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::Wall1LineType Discret::ELEMENTS::Wall1LineType::instance_;

Discret::ELEMENTS::Wall1LineType& Discret::ELEMENTS::Wall1LineType::instance() { return instance_; }

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 03/07|
  *----------------------------------------------------------------------*/
Discret::ELEMENTS::Wall1Line::Wall1Line(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::ELEMENTS::Wall1* parent, const int lline)
    : Core::Elements::FaceElement(id, owner)
{
  set_node_ids(nnode, nodeids);
  build_nodal_pointers(nodes);
  set_parent_master_element(parent, lline);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 03/07|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Wall1Line::Wall1Line(const Discret::ELEMENTS::Wall1Line& old)
    : Core::Elements::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Wall1Line::clone() const
{
  Discret::ELEMENTS::Wall1Line* newelement = new Discret::ELEMENTS::Wall1Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          farah 02/14 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Wall1Line::shape() const
{
  Core::FE::CellType distype_line = Core::FE::CellType::dis_none;

  switch (parent_master_element()->shape())
  {
    case Core::FE::CellType::tri3:
    {
      distype_line = Core::FE::CellType::line2;
      break;
    }
    case Core::FE::CellType::tri6:
    {
      distype_line = Core::FE::CellType::line3;
      break;
    }
    case Core::FE::CellType::quad4:
    {
      distype_line = Core::FE::CellType::line2;
      break;
    }
    case Core::FE::CellType::quad8:
    {
      distype_line = Core::FE::CellType::line3;
      break;
    }
    case Core::FE::CellType::quad9:
    {
      distype_line = Core::FE::CellType::line3;
      break;
    }
    case Core::FE::CellType::nurbs4:
    {
      distype_line = Core::FE::CellType::nurbs2;
      break;
    }
    case Core::FE::CellType::nurbs9:
    {
      distype_line = Core::FE::CellType::nurbs3;
      break;
    }
    default:
      FOUR_C_THROW("Discret::ELEMENTS::Wall1Line::Wall1Line: Unknown parent shape!");
  }

  return distype_line;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Wall1Line::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Wall1Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Wall1Line::unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this line element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              mgit 03/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Wall1Line::print(std::ostream& os) const
{
  os << "Wall1Line ";
  Element::print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
