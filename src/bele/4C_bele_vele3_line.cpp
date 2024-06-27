/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_bele_vele3.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN



Discret::ELEMENTS::Vele3LineType Discret::ELEMENTS::Vele3LineType::instance_;

Discret::ELEMENTS::Vele3LineType& Discret::ELEMENTS::Vele3LineType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Vele3Line::Vele3Line(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Core::Elements::Element* parent, const int lline)
    : Core::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lline);
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Vele3Line::Vele3Line(const Discret::ELEMENTS::Vele3Line& old)
    : Core::Elements::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Vele3Line::Clone() const
{
  Discret::ELEMENTS::Vele3Line* newelement = new Discret::ELEMENTS::Vele3Line(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Vele3Line::Shape() const
{
  switch (num_node())
  {
    case 2:
      return Core::FE::CellType::line2;
    case 3:
      return Core::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Vele3Line::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Vele3Line element does not support communication");

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Vele3Line::unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Vele3Line element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Vele3Line::print(std::ostream& os) const
{
  os << "Vele3Line ";
  Element::print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
