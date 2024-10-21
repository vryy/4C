// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_ele_calc_std.hpp"

#include "4C_fluid_ele_parameter_std.hpp"

FOUR_C_NAMESPACE_OPEN


template <Core::FE::CellType distype>
Discret::ELEMENTS::FluidEleCalcStd<distype>* Discret::ELEMENTS::FluidEleCalcStd<distype>::instance(
    Core::Utils::SingletonAction action)
{
  static auto singleton_owner = Core::Utils::make_singleton_owner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::FluidEleCalcStd<distype>>(
            new Discret::ELEMENTS::FluidEleCalcStd<distype>());
      });

  return singleton_owner.instance(action);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::FluidEleCalcStd<distype>::FluidEleCalcStd()
    : Discret::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc()
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

// template classes
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::tet10>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::wedge15>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::pyramid5>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::nurbs9>;
template class Discret::ELEMENTS::FluidEleCalcStd<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
