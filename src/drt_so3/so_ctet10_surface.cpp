/*!----------------------------------------------------------------------**
\file so_ctet10_surface.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
writen by : Alexander Volf
			alexander.volf@mytum.de  
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOTET
#ifdef CCADISCRET

#include "so_ctet10.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::Soctet10Surface::Soctet10Surface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::So_ctet10* parent,
                              const int lsurface) :
DRT::Element(id,element_soctet10surface,owner),
parent_(parent),
lsurface_(lsurface)
{

  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         maf 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soctet10Surface::Soctet10Surface(const DRT::ELEMENTS::Soctet10Surface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{

  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Soctet10Surface::Clone() const
{

  DRT::ELEMENTS::Soctet10Surface* newelement = new DRT::ELEMENTS::Soctet10Surface(*this);

  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Soctet10Surface::Shape() const
{
  return tri6;
}

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soctet10Surface::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Sote10Surface element does not support communication");

  return;
}

/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soctet10Surface::Unpack(const vector<char>& data)
{
  dserror("this Soctet10Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------***
 |  dtor (public)                                              maf 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soctet10Surface::~Soctet10Surface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soctet10Surface::Print(ostream& os) const
{
  os << "Soctet10Surface ";

  Element::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOTET
