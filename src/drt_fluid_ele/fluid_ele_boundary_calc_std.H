/*!----------------------------------------------------------------------
\file fluid_ele_boundary_calc_std.cpp

\brief evaluation of fluid terms at integration points

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
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

    template <DRT::Element::DiscretizationType distype>
    class FluidEleBoundaryCalcStd : public FluidBoundaryImpl<distype>
    {
      typedef DRT::ELEMENTS::FluidBoundaryImpl<distype> my;

     public:
      /// Singleton access method
      static FluidEleBoundaryCalcStd<distype>* Instance(bool create = true);

      /// Empty destructor
      virtual ~FluidEleBoundaryCalcStd(){};

      /// called upon destruction
      void Done();

     private:
      /// private constructor since we are singleton
      FluidEleBoundaryCalcStd();

    };  // class FluidEleBoundaryCalcStd

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
