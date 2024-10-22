// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
          Core::Utils::SingletonAction action = Core::Utils::SingletonAction::create);

     private:
      /// private constructor, since we are a Singleton.
      FluidEleCalcStd();
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
