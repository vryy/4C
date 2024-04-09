/*----------------------------------------------------------------------*/
/*! \file

\brief standard routines for calculation of fluid element

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_STD_HPP
#define FOUR_C_FLUID_ELE_CALC_STD_HPP

#include "baci_config.hpp"

#include "baci_fluid_ele_calc.hpp"
#include "baci_utils_singleton_owner.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype>
    class FluidEleCalcStd : public FluidEleCalc<distype>
    {
      typedef DRT::ELEMENTS::FluidEleCalc<distype> my;

     public:
      /// Singleton access method
      static FluidEleCalcStd<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

     private:
      /// private constructor, since we are a Singleton.
      FluidEleCalcStd();
    };
  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
