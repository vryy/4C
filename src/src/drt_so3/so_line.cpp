/*!----------------------------------------------------------------------
\file so_line.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_line.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                              gee 04/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::StructuralLine::StructuralLine(int id, int owner,
                                              int nnode, const int* nodeids,
                                              DRT::Node** nodes,
                                              DRT::Element* parent,
                                              const int lline) :
DRT::Element(id,element_structuralline,owner),
parent_(parent),
lline_(lline)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  // type of gaussian integration
  switch(Shape())
  {
  case line2:
    gaussrule_ = DRT::UTILS::intrule_line_2point;
  break;
  case line3:
    gaussrule_ = DRT::UTILS::intrule_line_3point;
  break;
  default:
      dserror("shape type unknown!\n");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 04/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::StructuralLine::StructuralLine(const DRT::ELEMENTS::StructuralLine& old) :
DRT::Element(old),
parent_(old.parent_),
lline_(old.lline_),
gaussrule_(old.gaussrule_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               gee 04/08|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::StructuralLine::Clone() const
{
  DRT::ELEMENTS::StructuralLine* newelement = new DRT::ELEMENTS::StructuralLine(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             gee 04/08|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::StructuralLine::Shape() const
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
 |  Pack data                                                  gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralLine::Pack(vector<char>& data) const
{
  dserror("StructuralLine element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralLine::Unpack(const vector<char>& data)
{
  dserror("StructuralLine element does not support communication");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                               gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralLine::Print(ostream& os) const
{
  os << "StructuralLine ";
  Element::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOLID3
