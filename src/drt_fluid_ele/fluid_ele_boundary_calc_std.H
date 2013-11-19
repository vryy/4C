/*!----------------------------------------------------------------------
\file fluid_ele_boundary_calc_Std.H
\brief

evaluate boundary conditions not requiring parent-element evaluations

<pre>
Maintainers: Anh-Tu Vuong & Andreas Rauch
             {vuong,rauch}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15264/-15240
</pre>
*----------------------------------------------------------------------*/
#ifndef FLUID_ELE_BOUNDARY_CALC_STD_H
#define FLUID_ELE_BOUNDARY_CALC_STD_H

#include "fluid_ele_boundary_calc.H"

namespace DRT
{
class Condition;
class Discretization;

namespace ELEMENTS
{
class FluidBoundary;

template<DRT::Element::DiscretizationType distype>
class FluidEleBoundaryCalcStd : public FluidBoundaryImpl<distype>
{
  typedef DRT::ELEMENTS::FluidBoundaryImpl<distype> my;

public:

  /// Singleton access method
  static FluidEleBoundaryCalcStd<distype> * Instance( bool create=true );

  /// Empty destructor
  virtual ~FluidEleBoundaryCalcStd() {};

  /// called upon destruction
  virtual void Done();

private:
  /// private constructor since we are singleton
  FluidEleBoundaryCalcStd();

}; // class FluidEleBoundaryCalcStd

} //namsepace ELEMENTS
} // namespace DRT

#endif
