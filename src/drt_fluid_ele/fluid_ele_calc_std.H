/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_std.H

\brief standard routines for calculation of fluid element

\level 1

<pre>
\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_STD_H
#define FLUID_ELE_CALC_STD_H

#include "fluid_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class FluidEleCalcStd : public FluidEleCalc<distype>
{
  typedef DRT::ELEMENTS::FluidEleCalc<distype> my;

public:
  virtual ~FluidEleCalcStd() {}

  /// Singleton access method
  static FluidEleCalcStd<distype>* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

private:
  /// private constructor, since we are a Singleton.
  FluidEleCalcStd();


};
}
}

#endif

