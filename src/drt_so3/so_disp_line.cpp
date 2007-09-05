/*!----------------------------------------------------------------------
\file fluid3_line.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger (Ursula)
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "so_disp.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispLine::SoDispLine(  int id,
                                        int owner,
                                        int nnode,
                                        const int* nodeids,
                                        DRT::Node** nodes,
                                        DRT::Elements::SoDispSurface* surfaceParent,
                                        DRT::Elements::SoDisp* parent,  
                                        const int lline) :
DRT::Element(id,element_sodispline,owner),
surfaceParent_(surfaceParent),
parent_(parent),
lline_(lline)
{
    SetNodeIds(nnode,nodeids);
    BuildNodalPointers(nodes);
    return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispLine::SoDispLine(const DRT::Elements::SoDispLine& old) :
DRT::Element(old),
surfaceParent_(old.surfaceParent_),
parent_(old.parent_),
lline_(old.lline_)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::SoDispLine::Clone() const
{
  DRT::Elements::SoDispLine* newelement = new DRT::Elements::SoDispLine(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::SoDispLine::Shape() const
{
  switch (NumNode())
  {
  case 2: return line2;
  case 3: return line3;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispLine::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this SoDispLine element does not support communication");

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispLine::Unpack(const vector<char>& data)
{
  dserror("this SoDispLine element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispLine::~SoDispLine()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispLine::Print(ostream& os) const
{
  os << "SoDispLine ";
  Element::Print(os);
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
