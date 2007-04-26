/*!----------------------------------------------------------------------
\file so_hex8_line.cpp
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
#ifdef TRILINOS_PACKAGE

#include "so_hex8.H"
#include "../discret/linalg_utils.H"
#include "../discret/drt_utils.H"
#include "../discret/drt_discret.H"
#include "../discret/drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../discret/dstrc.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Soh8Line::Soh8Line(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::So_hex8* parent,
                              const int lline) :
DRT::Element(id,element_soh8line,owner),
parent_(parent),
lline_(lline)
{
  DSTraceHelper dst("Soh8Line::Soh8Line");
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Soh8Line::Soh8Line(const DRT::Elements::Soh8Line& old) :
DRT::Element(old),
parent_(old.parent_),
lline_(old.lline_)
{
  DSTraceHelper dst("Soh8Line::Soh8Line");
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Soh8Line::Clone() const
{
  DSTraceHelper dst("Soh8Line::Clone");
  DRT::Elements::Soh8Line* newelement = new DRT::Elements::Soh8Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Soh8Line::Shape() const
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
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Soh8Line::Pack(vector<char>& data) const
{
  DSTraceHelper dst("Soh8Line::Pack");
  data.resize(0);

  dserror("this Soh8Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Soh8Line::Unpack(const vector<char>& data)
{
  DSTraceHelper dst("Soh8Line::Unpack");
  dserror("this line element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Soh8Line::~Soh8Line()
{
  DSTraceHelper dst("Soh8Line::~Soh8Line");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Soh8Line::Print(ostream& os) const
{
  DSTraceHelper dst("Soh8Line::Print");
  os << "Soh8Line ";
  Element::Print(os);
  return;
}

/*-----------------------------------------------------------------------*
 * Integrate a Line Neumann boundary condition (public)         maf 04/07*
 * ----------------------------------------------------------------------*/
int DRT::Elements::Soh8Line::EvaluateNeumann(ParameterList&         params,
                                             DRT::Discretization&   discretization,
                                             DRT::Condition&        condition,
                                             vector<int>&           lm,
                                             Epetra_SerialDenseVector& elevec1)
{
    return 0;
}                                             
                                            



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOH8
