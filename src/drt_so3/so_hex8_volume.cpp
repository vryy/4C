/*!----------------------------------------------------------------------
\file so_hex8_volume.cpp
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
DRT::Elements::Soh8Volume::Soh8Volume(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::So_hex8* parent,
                              const int lvolume) :
DRT::Element(id,element_soh8volume,owner),
parent_(parent),
lvolume_(lvolume)
{
  DSTraceHelper dst("Soh8Volume::Soh8Volume");
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Soh8Volume::Soh8Volume(const DRT::Elements::Soh8Volume& old) :
DRT::Element(old),
parent_(old.parent_),
lvolume_(old.lvolume_)
{
  DSTraceHelper dst("Soh8Volume::Soh8Volume");
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Soh8Volume::Clone() const
{
  DSTraceHelper dst("Soh8Volume::Clone");
  DRT::Elements::Soh8Volume* newelement = new DRT::Elements::Soh8Volume(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Soh8Volume::Shape() const
{
  return hex8;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Soh8Volume::Pack(vector<char>& data) const
{
  DSTraceHelper dst("Soh8Volume::Pack");
  data.resize(0);
  dserror("this Soh8Volume element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Soh8Volume::Unpack(const vector<char>& data)
{
  DSTraceHelper dst("Soh8Volume::Unpack");
  dserror("this Soh8Volume element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Soh8Volume::~Soh8Volume()
{
  DSTraceHelper dst("Soh8Volume::~Soh8Volume");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Soh8Volume::Print(ostream& os) const
{
  DSTraceHelper dst("Soh8Volume::Print");
  os << "Soh8Volume ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 * Integrate a Volume Neumann boundary condition (public)     maf 04/07*
 * ---------------------------------------------------------------------*/
int DRT::Elements::Soh8Volume::EvaluateNeumann(ParameterList&      params,
                                               Discretization&     discretization,
                                               Condition&          condition,
                                               vector<int>&        lm,
                                               Epetra_SerialDenseVector& elevec1)
{
    return 0;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOH8
