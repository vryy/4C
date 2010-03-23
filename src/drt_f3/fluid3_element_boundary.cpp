/*!----------------------------------------------------------------------
\file fluid3_surface.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3Boundary::Fluid3Boundary(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::Fluid3* parent,
                              const int lsurface) :
DRT::Element(id,element_fluid3boundary,owner),
parent_(parent),
lsurface_(lsurface)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3Boundary::Fluid3Boundary(const DRT::ELEMENTS::Fluid3Boundary& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Fluid3Boundary::Clone() const
{
  DRT::ELEMENTS::Fluid3Boundary* newelement = new DRT::ELEMENTS::Fluid3Boundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Fluid3Boundary::Shape() const
{
  switch (NumNode())
  {
  // 2D parent element -> Fluid3Boundary element: line2 and line3

  // Fluid3Boundary element: line2
  case 2:
    if(parent_->Shape() == quad4 || parent_->Shape() == tri3)
      return line2;
    else if (parent_->Shape() == nurbs4)
      return nurbs2;
    // 1D line element in a 3D volume
    else if(parent_->Shape() == hex8 ||
        parent_->Shape() == tet4 ||
        parent_->Shape() == wedge6 ||
        parent_->Shape() == pyramid5)
      return line2;
    // 1D line element in a 3D volume
    else if (parent_->Shape() == nurbs8)
      return nurbs2;
    else dserror("%d nodes of the Fluid3Boundary element does not fit to the distype %s of the parent element",
        NumNode(), DistypeToString(parent_->Shape()).c_str());

  // Fluid3Boundary element: line3
  case 3:
    if ((parent_->Shape() == quad8) || (parent_->Shape() == quad9))
      return line3;
    else if (parent_->Shape() == nurbs9)
      return nurbs3;
    // 1D line element in a 3D volume
    else if (parent_->Shape() == hex20 ||
        parent_->Shape() == hex27 ||
        parent_->Shape() == tet10 ||
        parent_->Shape() == wedge15)
      return line3;

  // Fluid3Boundary element: tri3 (surface)
    else if(parent_->Shape() == tet4 || parent_->Shape() == wedge6 || parent_->Shape() == pyramid5)
      return tri3;
    else dserror("%d nodes of the Fluid3Boundary element does not fit to the distype %s of the parent element",
        NumNode(), DistypeToString(parent_->Shape()).c_str());

  // Fluid3Boundary element: quad4
  case 4:
    if(parent_->Shape() == hex8 || parent_->Shape() == wedge6 || parent_->Shape() == pyramid5 )
      return quad4;
    else if (parent_->Shape() == nurbs8)
      return nurbs4;
    else dserror("%d nodes of the Fluid3Boundary element does not fit to the distype %s of the parent element",
        NumNode(), DistypeToString(parent_->Shape()).c_str());

  // Fluid3Boundary element: tri6
  case 6:
    if (parent_->Shape() == tet10 || parent_->Shape() == wedge15)
      return tri6;
    else dserror("%d nodes of the Fluid3Boundary element does not fit to the distype %s of the parent element",
        NumNode(), DistypeToString(parent_->Shape()).c_str());

  // Fluid3Boundary element: quad8
  case 8:
    if(parent_->Shape() == hex20 || parent_->Shape() == wedge15)
      return quad8;
    else dserror("%d nodes of the Fluid3Boundary element does not fit to the distype %s of the parent element",
        NumNode(), DistypeToString(parent_->Shape()).c_str());

  // Fluid3Boundary element: quad9
  case 9:
    if(parent_->Shape() == hex27)
        return quad9;
    else if (parent_->Shape() == nurbs27)
      return nurbs9;
    else dserror("%d nodes of the Fluid3Boundary element does not fit to the distype %s of the parent element",
        NumNode(), DistypeToString(parent_->Shape()).c_str());
  default:
    dserror("unexpected number of nodes %d for boundary element", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Boundary::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Fluid3Boundary element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Boundary::Unpack(const vector<char>& data)
{
  dserror("this Fluid3Boundary element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3Boundary::~Fluid3Boundary()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Boundary::Print(ostream& os) const
{
  os << "Fluid3Boundary ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Fluid3Boundary::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of Fluid3Boundary not implemented");
  vector<RCP<DRT::Element> > lines(0);
  return lines;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Fluid3Boundary::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of Fluid3Boundary not implemented");
  vector<RCP<DRT::Element> > surfaces(0);
  return surfaces;
}

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
