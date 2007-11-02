/*!----------------------------------------------------------------------
\file wall1_line.cpp
\brief

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15225
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET

#include "wall1.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 03/07|
  *----------------------------------------------------------------------*/
DRT::Elements::Wall1Line::Wall1Line(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::Wall1* parent,
                              const int lline) :
DRT::Element(id,element_wall1line,owner),
parent_(parent),
lline_(lline)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Line::Wall1Line(const DRT::Elements::Wall1Line& old) :
DRT::Element(old),
parent_(old.parent_),
lline_(old.lline_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Wall1Line::Clone() const
{
  DRT::Elements::Wall1Line* newelement = new DRT::Elements::Wall1Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Wall1Line::Shape() const
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
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1Line::Pack(vector<char>& data) const
{
  data.resize(0);
  
  dserror("this Wall1Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1Line::Unpack(const vector<char>& data)
{
  dserror("this line element does not support communication");
  return;
} 

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Wall1Line::~Wall1Line()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1Line::Print(ostream& os) const
{
  os << "Wall1Line ";
  Element::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_WALL1
