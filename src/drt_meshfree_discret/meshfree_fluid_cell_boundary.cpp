/*!---------------------------------------------------------------------------

\file meshfree_fluid_cell_boundary.cpp

\brief fluid cell boundary for meshfree discretisations

\maintainer Keijo Nissen

\level 3

*---------------------------------------------------------------------------*/

#include "meshfree_fluid_cell.H"
#include "drt_meshfree_node.H"
#include "../drt_lib/drt_utils_factory.H"

/*==========================================================================*\
 *                                                                          *
 * class MeshfreeFluidBoundaryType                                          *
 *                                                                          *
\*==========================================================================*/

/*--------------------------------------------------------------------------*
 |  create instance of MeshfreeFluidBoundaryType         (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeFluidBoundaryType DRT::ELEMENTS::MeshfreeFluidBoundaryType::instance_;

DRT::ELEMENTS::MeshfreeFluidBoundaryType& DRT::ELEMENTS::MeshfreeFluidBoundaryType::Instance()
{
  return instance_;
}

/*--------------------------------------------------------------------------*
 |  create object of MeshfreeFluidBoundaryType           (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::MeshfreeFluidBoundaryType::Create( const int id, const int owner )
{
  return Teuchos::null;
}


/*==========================================================================*\
 *                                                                          *
 * class MeshfreeFluidBoundary                                              *
 *                                                                          *
\*==========================================================================*/

/*--------------------------------------------------------------------------*
 |  ctor                                                 (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeFluidBoundary::MeshfreeFluidBoundary(
  int id,
  int owner,
  int npoint,
  int const * pointids,
  DRT::Node** points,
  MeshfreeFluid* parent,
  const int lsurface) :
  DRT::MESHFREE::Cell<DRT::FaceElement>(id,owner)
{
  SetPointIds(npoint,pointids);
  BuildPointPointers(points);
  SetParentMasterElement(parent,lsurface);

  // temporary assignement of nodes for call in DRT::Discretization::BuildLinesinCondition)
  // must and will be redefined in Face::AssignNodesToCells
  SetNodeIds(npoint,pointids);
  BuildNodalPointers(points);
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                            (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeFluidBoundary::MeshfreeFluidBoundary(const DRT::ELEMENTS::MeshfreeFluidBoundary& old) :
  DRT::MESHFREE::Cell<DRT::FaceElement>(old)
{
  return;
}

/*--------------------------------------------------------------------------*
 | returns pointer to deep copy of MeshfreeTransport     (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::MeshfreeFluidBoundary::Clone() const
{
  DRT::ELEMENTS::MeshfreeFluidBoundary* newelement = new DRT::ELEMENTS::MeshfreeFluidBoundary(*this);
  return newelement;
}

/*--------------------------------------------------------------------------*
 |  dtor                                                 (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeFluidBoundary::~MeshfreeFluidBoundary()
{
  return;
}

/*--------------------------------------------------------------------------*
 |  Return the shape of a MeshfreeFluidBoundary cell     (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::MeshfreeFluidBoundary::Shape() const
{
  return DRT::UTILS::getShapeOfBoundaryElement(NumPoint(), ParentElement()->Shape());
}

/*---------------------------------------------------------------------------*
 |  Return number of lines of boundary cell               (public) nis Feb12 |
 *---------------------------------------------------------------------------*/
inline int DRT::ELEMENTS::MeshfreeFluidBoundary::NumLine() const
{
  return DRT::UTILS::getNumberOfElementLines(Shape());
}

/*--------------------------------------------------------------------------*
 |  Return number of surfaces of boundary cell            (public) nis Feb12 |
 *---------------------------------------------------------------------------*/
inline int DRT::ELEMENTS::MeshfreeFluidBoundary::NumSurface() const
{
  return DRT::UTILS::getNumberOfElementSurfaces(Shape());
}

/*--------------------------------------------------------------------------*
 |  Get vector of Teuchos::RCPs to the lines (dummy)              (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::MeshfreeFluidBoundary::Lines()
{
  // surfaces, lines, and points have to be created by parent element
  dserror("Lines of MeshfreeTransportBoundary not implemented");

  // this is done to prevent compiler from moaning
  std::vector<Teuchos::RCP<DRT::Element> > lines(0);
  return lines;
}

/*--------------------------------------------------------------------------*
 |  Get vector of Teuchos::RCPs to the surfaces (dummy)           (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::MeshfreeFluidBoundary::Surfaces()
{
  // surfaces, lines, and points have to be created by parent element
  dserror("Surfaces of MeshfreeTransportBoundary not implemented");

  // this is done to prevent compiler from moaning
  std::vector<Teuchos::RCP<DRT::Element> > surfaces(0);
  return surfaces;
}

/*--------------------------------------------------------------------------*
 |  Pack data (dummy)                                    (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluidBoundary::Pack(DRT::PackBuffer& data) const
{
  dserror("This MeshfreeFluidBoundary cell does not support communication.");

  return;
}

/*--------------------------------------------------------------------------*
 |  Unpack data (dummy)                                  (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluidBoundary::Unpack(const std::vector<char>& data)
{
  dserror("This MeshfreeFluidBoundary cell doesnot support communication.");
  return;
}


/*--------------------------------------------------------------------------*
 |  Print this cell                                      (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluidBoundary::Print(std::ostream& os) const
{
  os << "MeshfreeFluidBoundary ";
  Element::Print(os);
  return;
}
