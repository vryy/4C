/*!----------------------------------------------------------------------**
\file so_ctet10_line.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOTET
#ifdef CCADISCRET

#include "so_ctet10.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Soctet10Line::Soctet10Line(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::So_ctet10* parent,
                              const int lline) :
DRT::Element(id,element_soctet10line,owner),
parent_(parent),
lline_(lline)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Soctet10Line::Soctet10Line(const DRT::Elements::Soctet10Line& old) :
DRT::Element(old),
parent_(old.parent_),
lline_(old.lline_)
{
  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Soctet10Line::Clone() const
{

  DRT::Elements::Soctet10Line* newelement = new DRT::Elements::Soctet10Line(*this);
  return newelement;
}

/*----------------------------------------------------------------------***
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Soctet10Line::Shape() const
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
void DRT::Elements::Soctet10Line::Pack(vector<char>& data) const
{
  data.resize(0);

  dserror("this Soctet10Line element does not support communication");

  return;
}

/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Soctet10Line::Unpack(const vector<char>& data)
{
  dserror("this line element does not support communication");
  return;
}

/*----------------------------------------------------------------------***
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Soctet10Line::~Soctet10Line()
{
  return;
}


/*----------------------------------------------------------------------***
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Soctet10Line::Print(ostream& os) const
{
  os << "Soctet10Line ";

  Element::Print(os);
  return;
}

/*-----------------------------------------------------------------------***
 * Integrate a Line Neumann boundary condition (public)         maf 04/07*
 * ----------------------------------------------------------------------*/
int DRT::Elements::Soctet10Line::EvaluateNeumann(ParameterList&         params,
                                             DRT::Discretization&   discretization,
                                             DRT::Condition&        condition,
                                             vector<int>&           lm,
                                             Epetra_SerialDenseVector& elevec1)
{
  dserror("Neumann condition on line not implemented");
   
  return 0;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOTET
