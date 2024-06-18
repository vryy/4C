/*----------------------------------------------------------------------*/
/*! \file

\brief dummy 3D boundary element without any physics


\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_bele_bele3.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::Bele3LineType Discret::ELEMENTS::Bele3LineType::instance_;

Discret::ELEMENTS::Bele3LineType& Discret::ELEMENTS::Bele3LineType::Instance() { return instance_; }


/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Bele3Line::Bele3Line(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::ELEMENTS::Bele3* parent, const int lline)
    : Core::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lline);
  set_num_dof_per_node(parent->NumDofPerNode(*nodes[0]));
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Bele3Line::Bele3Line(const Discret::ELEMENTS::Bele3Line& old)
    : Core::Elements::FaceElement(old), numdofpernode_(old.numdofpernode_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Bele3Line::Clone() const
{
  Discret::ELEMENTS::Bele3Line* newelement = new Discret::ELEMENTS::Bele3Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Bele3Line::Shape() const
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
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Bele3Line::pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Bele3Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Bele3Line::unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Bele3Line element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Bele3Line::Print(std::ostream& os) const
{
  os << "Bele3_" << numdofpernode_ << "Line ";
  Element::Print(os);
  return;
}

FOUR_C_NAMESPACE_CLOSE
