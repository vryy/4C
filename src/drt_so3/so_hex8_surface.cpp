/*!----------------------------------------------------------------------
\file so_hex8_surface.cpp
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

#include "so_hex8.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soh8Surface::Soh8Surface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::So_hex8* parent,
                              const int lsurface) :
DRT::Element(id,element_soh8surface,owner),
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
DRT::ELEMENTS::Soh8Surface::Soh8Surface(const DRT::ELEMENTS::Soh8Surface& old) :
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
DRT::Element* DRT::ELEMENTS::Soh8Surface::Clone() const
{
  DRT::ELEMENTS::Soh8Surface* newelement = new DRT::ELEMENTS::Soh8Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Soh8Surface::Shape() const
{
  return quad4;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Surface::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Soh8Surface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Surface::Unpack(const vector<char>& data)
{
  dserror("this Soh8Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Soh8Surface::~Soh8Surface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Soh8Surface::Print(ostream& os) const
{
  os << "Soh8Surface ";
  Element::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOH8
