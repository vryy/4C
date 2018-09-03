/*----------------------------------------------------------------------*/
/*!
\file vele3_surface.cpp

\brief volume element

\maintainer Jonas Eichinger

\level 2
*/
/*----------------------------------------------------------------------*/

#include "vele3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_factory.H"


DRT::ELEMENTS::Vele3SurfaceType DRT::ELEMENTS::Vele3SurfaceType::instance_;

DRT::ELEMENTS::Vele3SurfaceType& DRT::ELEMENTS::Vele3SurfaceType::Instance() { return instance_; }


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 05/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Surface::Vele3Surface(int id, int owner, int nnode, const int* nodeids,
    DRT::Node** nodes, DRT::ELEMENTS::Vele3* parent, const int lsurface)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lsurface);
  return;
}



/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Surface::Vele3Surface(const DRT::ELEMENTS::Vele3Surface& old)
    : DRT::FaceElement(old)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Vele3Surface::Clone() const
{
  DRT::ELEMENTS::Vele3Surface* newelement = new DRT::ELEMENTS::Vele3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Vele3Surface::Shape() const
{
  switch (NumNode())
  {
    case 3:
      return tri3;
    case 4:
      return quad4;
    case 6:
      return tri6;
    case 8:
      return quad8;
    case 9:
      return quad9;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Surface::Pack(DRT::PackBuffer& data) const
{
  dserror("this Vele3Surface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Surface::Unpack(const std::vector<char>& data)
{
  dserror("this Vele3Surface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Vele3Surface::~Vele3Surface() { return; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Vele3Surface::Print(std::ostream& os) const
{
  os << "Vele3Surface " << DRT::DistypeToString(Shape());
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Vele3Surface::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Vele3Line, Vele3Surface>(DRT::UTILS::buildLines, this);
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Vele3Surface::Surfaces()
{
  std::vector<Teuchos::RCP<DRT::Element>> surfaces(1);
  surfaces[0] = Teuchos::rcp(this, false);
  return surfaces;
}



/*----------------------------------------------------------------------*
 |  get optimal gauss rule                                   gammi 04/07|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule2D DRT::ELEMENTS::Vele3Surface::getOptimalGaussrule(
    const DRT::Element::DiscretizationType& distype) const
{
  DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule2D_undefined;
  switch (distype)
  {
    case DRT::Element::quad4:
      rule = DRT::UTILS::intrule_quad_4point;
      break;
    case DRT::Element::quad8:
    case DRT::Element::quad9:
      rule = DRT::UTILS::intrule_quad_9point;
      break;
    case DRT::Element::tri3:
      rule = DRT::UTILS::intrule_tri_3point;
      break;
    case DRT::Element::tri6:
      rule = DRT::UTILS::intrule_tri_6point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}
