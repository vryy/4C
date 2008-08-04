/*!----------------------------------------------------------------------
\file condif3_surface.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "condif3.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 06/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif3Surface::Condif3Surface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::Condif3* parent,
                              const int lsurface) :
DRT::Element(id,element_condif3surface,owner),
parent_(parent),
lsurface_(lsurface)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 06/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif3Surface::Condif3Surface(const DRT::ELEMENTS::Condif3Surface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it     (public) gjb 06/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Condif3Surface::Clone() const
{
  DRT::ELEMENTS::Condif3Surface* newelement = new DRT::ELEMENTS::Condif3Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                    (public)  gjb 06/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Condif3Surface::Shape() const
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
 |  Pack data (public)                                        gjb 06/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Surface::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this Condif3Surface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                      gjb 06/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Surface::Unpack(const vector<char>& data)
{
  dserror("this Condif3Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                             gjb 06/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif3Surface::~Condif3Surface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 06/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Surface::Print(ostream& os) const
{
  os << "Condif3Surface ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 06/08 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Condif3Surface::Lines()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Condif3Lines not implemented");
  vector<RCP<DRT::Element> > lines(0);
  return lines;
}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
