/*!----------------------------------------------------------------------**
\file so_tet4_surface.cpp
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
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_tet4.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

//extern "C"
//{
//#include "../headers/standardtypes.h"
//}
//#include "../drt_lib/dstrc.H"
 


/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/

DRT::ELEMENTS::Sotet4Surface::Sotet4Surface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::So_tet4* parent,
                              const int lsurface) :
DRT::Element(id,element_sotet4surface,owner),
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
DRT::ELEMENTS::Sotet4Surface::Sotet4Surface(const DRT::ELEMENTS::Sotet4Surface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{
  //DSTraceHelper dst("Sotet4Surface::Sotet4Surface");
  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Sotet4Surface::Clone() const
{
  //DSTraceHelper dst("Sotet4Surface::Clone");
  DRT::ELEMENTS::Sotet4Surface* newelement = new DRT::ELEMENTS::Sotet4Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Sotet4Surface::Shape() const
{

  return tri6;
}

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sotet4Surface::Pack(vector<char>& data) const
{
  //DSTraceHelper dst("Sotet4Surface::Pack");
  data.resize(0);
  dserror("this Sote10Surface element does not support communication");

  return;
}

/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sotet4Surface::Unpack(const vector<char>& data)
{
  //DSTraceHelper dst("Sotet4Surface::Unpack");
 
  dserror("this Sotet4Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------***
 |  dtor (public)                                              maf 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sotet4Surface::~Sotet4Surface()
{
  //DSTraceHelper dst("Sotet4Surface::~Sotet4Surface");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sotet4Surface::Print(ostream& os) const
{
  //DSTraceHelper dst("Sotet4Surface::Print");
  os << "Sotet4Surface ";
  Element::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOLID3
