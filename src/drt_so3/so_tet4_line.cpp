/*!----------------------------------------------------------------------**
\file so_tet4_line.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOTET4
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "so_tet4.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../drt_lib/dstrc.H"



/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Sotet4Line::Sotet4Line(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::So_tet4* parent,
                              const int lline) :
DRT::Element(id,element_sotet4line,owner),
parent_(parent),
lline_(lline)
{
  DSTraceHelper dst("Sotet4Line::Sotet4Line");
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Sotet4Line::Sotet4Line(const DRT::Elements::Sotet4Line& old) :
DRT::Element(old),
parent_(old.parent_),
lline_(old.lline_)
{
  DSTraceHelper dst("Sotet4Line::Sotet4Line");
  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Sotet4Line::Clone() const
{
  DSTraceHelper dst("Sotet4Line::Clone");
  DRT::Elements::Sotet4Line* newelement = new DRT::Elements::Sotet4Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------***
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Sotet4Line::Shape() const
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

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Sotet4Line::Pack(vector<char>& data) const
{
  DSTraceHelper dst("Sotet4Line::Pack");
  data.resize(0);

  dserror("this Sotet4Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Sotet4Line::Unpack(const vector<char>& data)
{
  DSTraceHelper dst("Sotet4Line::Unpack");
  dserror("this line element does not support communication");
  return;
}

/*----------------------------------------------------------------------***
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Sotet4Line::~Sotet4Line()
{
  DSTraceHelper dst("Sotet4Line::~Sotet4Line");
  return;
}


/*----------------------------------------------------------------------***
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Sotet4Line::Print(ostream& os) const
{
  DSTraceHelper dst("Sotet4Line::Print");
  os << "Sote10Line ";
  Element::Print(os);
  return;
}

/*-----------------------------------------------------------------------***
 * Integrate a Line Neumann boundary condition (public)         maf 04/07*
 * ----------------------------------------------------------------------*/
int DRT::Elements::Sotet4Line::EvaluateNeumann(ParameterList&         params,
                                             DRT::Discretization&   discretization,
                                             DRT::Condition&        condition,
                                             vector<int>&           lm,
                                             Epetra_SerialDenseVector& elevec1)
{
  dserror("Neumann condition on line not implemented");
  return 0;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOTET4
