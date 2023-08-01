/*----------------------------------------------------------------------*/
/*! \file

\brief standard routines for calculation of fluid element

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_CALC_STD_H
#define BACI_FLUID_ELE_CALC_STD_H

#include "baci_fluid_ele_calc.H"
#include "baci_utils_singleton_owner.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype>
    class FluidEleCalcStd : public FluidEleCalc<distype>
    {
      typedef DRT::ELEMENTS::FluidEleCalc<distype> my;

     public:
      virtual ~FluidEleCalcStd() {}

      /// Singleton access method
      static FluidEleCalcStd<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

     private:
      /// private constructor, since we are a Singleton.
      FluidEleCalcStd();
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
