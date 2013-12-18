/*----------------------------------------------------------------------*/
/*!
\file meshfree_fluid_cell_calc_std.H

\brief standard routines for calculation of meshfree fluid cell

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MESHFREE_FLUID_CELL_CALC_STD_H
#define MESHFREE_FLUID_CELL_CALC_STD_H

#include "meshfree_fluid_cell_calc.H"

namespace DRT
{
namespace ELEMENTS
{
  template<DRT::Element::DiscretizationType distype>
  class MeshfreeFluidCellCalcStd : public MeshfreeFluidCellCalc<distype>
  {
  private:
    /// private constructor, since we are a Singleton.
    MeshfreeFluidCellCalcStd();

    typedef MeshfreeFluidCellCalcStd<distype> my;

  public:

    /// destructor
    virtual ~MeshfreeFluidCellCalcStd() {}

    /// Singleton access method
    static MeshfreeFluidCellCalcStd<distype> * Instance( bool create=true );

    /// called upon destruction
    virtual void Done();
  };
}
}

#endif

