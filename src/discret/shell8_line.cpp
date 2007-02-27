/*!----------------------------------------------------------------------
\file drt_elementsurface.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "shell8.H"
#include "linalg_utils.H"
#include "drt_utils.H"
#include "drt_discret.H"
#include "drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
#include "../shell8/shell8.h"
}
#include "dstrc.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Line::Shell8Line(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::Shell8* parent,
                              const int lline) :
DRT::Element(id,element_shell8line,owner),
parent_(parent),
lline_(lline)
{
  DSTraceHelper dst("Shell8Line::Shell8Line");
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Line::Shell8Line(const DRT::Elements::Shell8Line& old) :
DRT::Element(old),
parent_(old.parent_),
lline_(old.lline_)
{
  DSTraceHelper dst("Shell8Line::Shell8Line");
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Shell8Line::Clone() const
{
  DSTraceHelper dst("Shell8Line::Clone");
  DRT::Elements::Shell8Line* newelement = new DRT::Elements::Shell8Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8Line::Pack(vector<char>& data) const
{
  DSTraceHelper dst("Shell8Line::Pack");  
  data.resize(0);
  
  dserror("this Shell8Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8Line::Unpack(const vector<char>& data)
{
  DSTraceHelper dst("Shell8Line::Unpack");  
  dserror("this line element does not support communication");
  return;
} 

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Line::~Shell8Line()
{
  DSTraceHelper dst("Shell8Line::~Shell8Line");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8Line::Print(ostream& os) const
{
  DSTraceHelper dst("Shell8Line::Print");
  os << "Shell8Line ";
  Element::Print(os);
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SHELL8
