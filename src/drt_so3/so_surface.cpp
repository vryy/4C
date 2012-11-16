/*!----------------------------------------------------------------------
\file so_surface.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "so_surface.H"
#include "so_line.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"

DRT::ELEMENTS::StructuralSurfaceType DRT::ELEMENTS::StructuralSurfaceType::instance_;


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::StructuralSurfaceType::Create( const int id, const int owner )
{
  //return Teuchos::rcp( new StructuralSurface( id, owner ) );
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                              gee 04/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::StructuralSurface::StructuralSurface(int id, int owner,
                                                    int nnode, const int* nodeids,
                                                    DRT::Node** nodes,
                                                    DRT::Element* parent,
                                                    const int lsurface) :
DRT::Element(id,owner),
parent_(parent),
lsurface_(lsurface),
gaussrule_(DRT::UTILS::intrule2D_undefined)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  // type of gaussian integration
  switch(Shape())
  {
  case tri3:
    gaussrule_ = DRT::UTILS::intrule_tri_3point;
  break;
  case tri6:
    gaussrule_ = DRT::UTILS::intrule_tri_6point;
  break;
  case quad4:
    gaussrule_ = DRT::UTILS::intrule_quad_4point;
  break;
  case quad8:
    gaussrule_ = DRT::UTILS::intrule_quad_9point;
  break;
  case quad9:
    gaussrule_ = DRT::UTILS::intrule_quad_9point;
  break;
  case nurbs9:
    gaussrule_ = DRT::UTILS::intrule_quad_9point;
  break;
  default:
      dserror("shape type unknown!\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 04/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::StructuralSurface::StructuralSurface(const DRT::ELEMENTS::StructuralSurface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_),
gaussrule_(old.gaussrule_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 04/08|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::StructuralSurface::Clone() const
{
  DRT::ELEMENTS::StructuralSurface* newelement = new DRT::ELEMENTS::StructuralSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::StructuralSurface::Shape() const
{
  switch (NumNode())
  {
  case 3: return tri3;
  case 6: return tri6;
  case 4: return quad4;
  case 8: return quad8;
  case 9: return quad9;
  default: dserror("Unknown shape of surface element (unknown number of nodes)");
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::Pack(DRT::PackBuffer& data) const
{
  dserror("this StructuralSurface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::Unpack(const std::vector<char>& data)
{
  dserror("this StructuralSurface element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::Print(ostream& os) const
{
  os << "StructuralSurface ";
  Element::Print(os);
  return;
}

std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::StructuralSurface::Lines()
{
   return DRT::UTILS::ElementBoundaryFactory
     <DRT::ELEMENTS::StructuralLine,DRT::ELEMENTS::StructuralSurface>(DRT::UTILS::buildLines,this);
}



