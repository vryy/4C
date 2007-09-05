/*!----------------------------------------------------------------------
\file so_disp_surface.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "so_disp.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"

//extern "C"
//{
//#include "../headers/standardtypes.h"
//}
//#include "../drt_lib/dstrc.H"
 


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispSurface::SoDispSurface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::SoDisp* parent,
                              const int lsurface) :
DRT::Element(id,element_sodispsurface,owner),
parent_(parent),
lsurface_(lsurface)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispSurface::SoDispSurface(const DRT::Elements::SoDispSurface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::SoDispSurface::Clone() const
{
  DRT::Elements::SoDispSurface* newelement = new DRT::Elements::SoDispSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::SoDispSurface::Shape() const
{
    switch (NumNode())
    {
    case 3: return tri3;
    case 4: return quad4;
    case 6: return tri6;
    case 8: return quad8;
    case 9: return quad9;
    default:
      dserror("unexpected number of nodes %d", NumNode());
    }
    return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispSurface::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this SoDispSurface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispSurface::Unpack(const vector<char>& data)
{
  dserror("this SoDispSurface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispSurface::~SoDispSurface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispSurface::Print(ostream& os) const
{
  os << "SoDispSurface ";
  Element::Print(os);
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOH8
