// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solid_poro_3D_ele_factory.hpp"

FOUR_C_NAMESPACE_OPEN



Discret::ELEMENTS::SolidPoroPressureBasedCalcVariant
Discret::ELEMENTS::create_solid_poro_pressure_based_calculation_interface(
    Core::FE::CellType celltype)
{
  return Core::FE::cell_type_switch<Discret::ELEMENTS::Internal::ImplementedSolidPoroCellTypes>(
      celltype, [&](auto celltype_t)
      { return create_solid_poro_pressure_based_calculation_interface<celltype_t()>(); });
}

template <Core::FE::CellType celltype>
Discret::ELEMENTS::SolidPoroPressureBasedCalcVariant
Discret::ELEMENTS::create_solid_poro_pressure_based_calculation_interface()
{
  return SolidPoroPressureBasedEleCalc<celltype>();
}

Discret::ELEMENTS::SolidPoroPressureVelocityBasedCalcVariant
Discret::ELEMENTS::create_solid_poro_pressure_velocity_based_calculation_interface(
    Core::FE::CellType celltype)
{
  return Core::FE::cell_type_switch<
      Discret::ELEMENTS::Internal::ImplementedSolidPoroPressureVelocityBasedCellTypes>(celltype,
      [&](auto celltype_t)
      { return create_solid_poro_pressure_velocity_based_calculation_interface<celltype_t()>(); });
}


template <Core::FE::CellType celltype>
Discret::ELEMENTS::SolidPoroPressureVelocityBasedCalcVariant
Discret::ELEMENTS::create_solid_poro_pressure_velocity_based_calculation_interface()
{
  return SolidPoroPressureVelocityBasedEleCalc<celltype>();
}



FOUR_C_NAMESPACE_CLOSE
