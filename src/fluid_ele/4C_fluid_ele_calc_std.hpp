/*----------------------------------------------------------------------*/
/*! \file

\brief standard routines for calculation of fluid element

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_STD_HPP
#define FOUR_C_FLUID_ELE_CALC_STD_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_calc.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    template <Core::FE::CellType distype>
    class FluidEleCalcStd : public FluidEleCalc<distype>
    {
      typedef Discret::ELEMENTS::FluidEleCalc<distype> my;

     public:
      /// Singleton access method
      static FluidEleCalcStd<distype>* instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

     private:
      /// private constructor, since we are a Singleton.
      FluidEleCalcStd();
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
