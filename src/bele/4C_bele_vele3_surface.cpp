/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_bele_vele3.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_lib_discret.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::Vele3SurfaceType DRT::ELEMENTS::Vele3SurfaceType::instance_;

DRT::ELEMENTS::Vele3SurfaceType& DRT::ELEMENTS::Vele3SurfaceType::Instance() { return instance_; }


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 05/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Surface::Vele3Surface(int id, int owner, int nnode, const int* nodeids,
    CORE::Nodes::Node** nodes, DRT::ELEMENTS::Vele3* parent, const int lsurface)
    : CORE::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}



/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Surface::Vele3Surface(const DRT::ELEMENTS::Vele3Surface& old)
    : CORE::Elements::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::Vele3Surface::Clone() const
{
  DRT::ELEMENTS::Vele3Surface* newelement = new DRT::ELEMENTS::Vele3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Vele3Surface::Shape() const
{
  switch (num_node())
  {
    case 3:
      return CORE::FE::CellType::tri3;
    case 4:
      return CORE::FE::CellType::quad4;
    case 6:
      return CORE::FE::CellType::tri6;
    case 8:
      return CORE::FE::CellType::quad8;
    case 9:
      return CORE::FE::CellType::quad9;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Surface::Pack(CORE::COMM::PackBuffer& data) const
{
  FOUR_C_THROW("this Vele3Surface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Surface::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Vele3Surface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Surface::Print(std::ostream& os) const
{
  os << "Vele3Surface " << CORE::FE::CellTypeToString(Shape());
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::Vele3Surface::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<Vele3Line, Vele3Surface>(CORE::COMM::buildLines, *this);
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::Vele3Surface::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}



/*----------------------------------------------------------------------*
 |  get optimal gauss rule                                   gammi 04/07|
 *----------------------------------------------------------------------*/
CORE::FE::GaussRule2D DRT::ELEMENTS::Vele3Surface::get_optimal_gaussrule(
    const CORE::FE::CellType& distype) const
{
  CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::undefined;
  switch (distype)
  {
    case CORE::FE::CellType::quad4:
      rule = CORE::FE::GaussRule2D::quad_4point;
      break;
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
      rule = CORE::FE::GaussRule2D::quad_9point;
      break;
    case CORE::FE::CellType::tri3:
      rule = CORE::FE::GaussRule2D::tri_3point;
      break;
    case CORE::FE::CellType::tri6:
      rule = CORE::FE::GaussRule2D::tri_6point;
      break;
    default:
      FOUR_C_THROW("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

FOUR_C_NAMESPACE_CLOSE
