/*!----------------------------------------------------------------------
\file ale2_line.cpp
\brief

<pre>
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#ifdef CCADISCRET

#include "ale2.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Ale2Line::Ale2Line(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::Ale2* parent,
                              const int lline) :
DRT::Element(id,element_ale2line,owner),
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
DRT::Elements::Ale2Line::Ale2Line(const DRT::Elements::Ale2Line& old) :
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
DRT::Element* DRT::Elements::Ale2Line::Clone() const
{
  DRT::Elements::Ale2Line* newelement = new DRT::Elements::Ale2Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Ale2Line::Shape() const
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
void DRT::Elements::Ale2Line::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Ale2Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Ale2Line::Unpack(const vector<char>& data)
{
  dserror("this Ale2Line element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Ale2Line::~Ale2Line()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Ale2Line::Print(ostream& os) const
{
  os << "Ale2Line ";
  Element::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_ALE2
