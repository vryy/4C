/*!----------------------------------------------------------------------
\file fluid2_line.cpp
\brief

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "fluid2.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"




/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Line::Fluid2Line(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::Fluid2* parent,
                              const int lline) :
DRT::Element(id,element_fluid2line,owner),
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
DRT::ELEMENTS::Fluid2Line::Fluid2Line(const DRT::ELEMENTS::Fluid2Line& old) :
DRT::Element(old),
parent_(old.parent_),
lline_(old.lline_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Fluid2Line::Clone() const
{
  DRT::ELEMENTS::Fluid2Line* newelement = new DRT::ELEMENTS::Fluid2Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Fluid2Line::Shape() const
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
void DRT::ELEMENTS::Fluid2Line::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Fluid2Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Line::Unpack(const vector<char>& data)
{
  dserror("this Fluid2Line element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Line::~Fluid2Line()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Line::Print(ostream& os) const
{
  os << "Fluid2Line ";
  Element::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID2
