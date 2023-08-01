/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of fluid terms at integration points of boundaries

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_BOUNDARY_CALC_STD_H
#define BACI_FLUID_ELE_BOUNDARY_CALC_STD_H

#include "baci_fluid_ele_boundary_calc.H"

#include "baci_utils_singleton_owner.H"

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
      static FluidEleBoundaryCalcStd<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

     private:
      /// private constructor since we are singleton
      FluidEleBoundaryCalcStd();

    };  // class FluidEleBoundaryCalcStd

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
