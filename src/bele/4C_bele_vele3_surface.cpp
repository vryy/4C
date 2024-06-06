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


Discret::ELEMENTS::Vele3SurfaceType Discret::ELEMENTS::Vele3SurfaceType::instance_;

Discret::ELEMENTS::Vele3SurfaceType& Discret::ELEMENTS::Vele3SurfaceType::Instance()
{
  return instance_;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 05/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Vele3Surface::Vele3Surface(int id, int owner, int nnode, const int* nodeids,
    Core::Nodes::Node** nodes, Discret::ELEMENTS::Vele3* parent, const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}



/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Vele3Surface::Vele3Surface(const Discret::ELEMENTS::Vele3Surface& old)
    : Core::Elements::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Vele3Surface::Clone() const
{
  Discret::ELEMENTS::Vele3Surface* newelement = new Discret::ELEMENTS::Vele3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Vele3Surface::Shape() const
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
void Discret::ELEMENTS::Vele3Surface::Pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("this Vele3Surface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Vele3Surface::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("this Vele3Surface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Vele3Surface::Print(std::ostream& os) const
{
  os << "Vele3Surface " << Core::FE::CellTypeToString(Shape());
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Vele3Surface::Lines()
{
  return Core::Communication::ElementBoundaryFactory<Vele3Line, Vele3Surface>(
      Core::Communication::buildLines, *this);
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Vele3Surface::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}



/*----------------------------------------------------------------------*
 |  get optimal gauss rule                                   gammi 04/07|
 *----------------------------------------------------------------------*/
Core::FE::GaussRule2D Discret::ELEMENTS::Vele3Surface::get_optimal_gaussrule(
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
