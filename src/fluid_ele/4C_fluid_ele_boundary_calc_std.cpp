// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_ele_boundary_calc_std.hpp"

#include "4C_fem_general_elementtype.hpp"
#include "4C_fluid_ele_parameter_std.hpp"

FOUR_C_NAMESPACE_OPEN

template <Core::FE::CellType distype>
Discret::Elements::FluidEleBoundaryCalcStd<distype>*
Discret::Elements::FluidEleBoundaryCalcStd<distype>::instance(Core::Utils::SingletonAction action)
{
  static auto singleton_owner = Core::Utils::make_singleton_owner(
      []()
      {
        return std::unique_ptr<Discret::Elements::FluidEleBoundaryCalcStd<distype>>(
            new Discret::Elements::FluidEleBoundaryCalcStd<distype>());
      });

  return singleton_owner.instance(action);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::Elements::FluidEleBoundaryCalcStd<distype>::FluidEleBoundaryCalcStd()
    : Discret::Elements::FluidBoundaryImpl<distype>::FluidBoundaryImpl()
{
  // pointer to class FluidImplParameter
  my::fldpara_ = Discret::Elements::FluidEleParameterStd::instance();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::quad4>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::quad8>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::quad9>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::tri3>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::tri6>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::line2>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::line3>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::nurbs2>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::nurbs3>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::nurbs4>;
template class Discret::Elements::FluidEleBoundaryCalcStd<Core::FE::CellType::nurbs9>;

FOUR_C_NAMESPACE_CLOSE
