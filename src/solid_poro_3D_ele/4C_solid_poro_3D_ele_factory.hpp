/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_FACTORY_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_solid_3D_ele_factory_lib.hpp"
#include "4C_solid_poro_3D_ele_calc_pressure_based.hpp"
#include "4C_solid_poro_3D_ele_calc_pressure_velocity_based.hpp"

#include <memory>
#include <variant>

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{

  namespace Details
  {
    using ImplementedSolidPoroCellTypes = Core::FE::CelltypeSequence<Core::FE::CellType::hex8,
        Core::FE::CellType::hex27, Core::FE::CellType::tet4, Core::FE::CellType::tet10>;

    using PoroPressureBasedEvaluators =
        Core::FE::apply_celltype_sequence<Discret::ELEMENTS::SolidPoroPressureBasedEleCalc,
            ImplementedSolidPoroCellTypes>;

    using SolidPoroPressureBasedEvaluators = Core::FE::Join<PoroPressureBasedEvaluators>;

    using ImplementedSolidPoroPressureVelocityBasedCellTypes =
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8, Core::FE::CellType::hex27,
            Core::FE::CellType::tet4, Core::FE::CellType::tet10>;


    using PoroPressureVelocityBasedEvaluators =
        Core::FE::apply_celltype_sequence<Discret::ELEMENTS::SolidPoroPressureVelocityBasedEleCalc,
            ImplementedSolidPoroPressureVelocityBasedCellTypes>;


    using SolidPoroPressureVelocityBasedEvaluators =
        Core::FE::Join<PoroPressureVelocityBasedEvaluators>;


  }  // namespace Details


  using SolidPoroPressureBasedCalcVariant =
      CreateVariantType<Details::SolidPoroPressureBasedEvaluators>;

  SolidPoroPressureBasedCalcVariant create_solid_poro_pressure_based_calculation_interface(
      Core::FE::CellType celltype);

  template <Core::FE::CellType celltype>
  SolidPoroPressureBasedCalcVariant create_solid_poro_pressure_based_calculation_interface();

  using SolidPoroPressureVelocityBasedCalcVariant =
      CreateVariantType<Details::SolidPoroPressureVelocityBasedEvaluators>;


  SolidPoroPressureVelocityBasedCalcVariant
  create_solid_poro_pressure_velocity_based_calculation_interface(Core::FE::CellType celltype);

  template <Core::FE::CellType celltype>
  SolidPoroPressureVelocityBasedCalcVariant
  create_solid_poro_pressure_velocity_based_calculation_interface();

}  // namespace Discret::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif