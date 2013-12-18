/*----------------------------------------------------------------------*/
/*!
\file meshfree_fluid_cell_boundary_calc_std.H

\brief standard routines for calculation of meshfree fluid boundary cell

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MESHFREE_FLUID_CELL_BOUNDARY_CALC_STD_H
#define MESHFREE_FLUID_CELL_BOUNDARY_CALC_STD_H

#include "meshfree_fluid_cell_boundary_calc.H"

namespace DRT
{
class Condition;
class Discretization;

namespace ELEMENTS
{
class MeshfreeFluidBoundary;

template<DRT::Element::DiscretizationType distype>
class MeshfreeFluidBoundaryCalcStd : public MeshfreeFluidBoundaryCalc<distype>
{
  typedef DRT::ELEMENTS::MeshfreeFluidBoundaryCalc<distype> my;

public:

  /// Singleton access method
  static MeshfreeFluidBoundaryCalcStd<distype>* Instance( bool create=true );

  /// Empty destructor
  virtual ~MeshfreeFluidBoundaryCalcStd() {};

  /// called upon destruction
  void Done();

private:
  /// private constructor since we are singleton
  MeshfreeFluidBoundaryCalcStd();

}; // class MeshfreeFluidCellBoundaryCalcStd

} // namsepace ELEMENTS
} // namespace DRT

#endif
