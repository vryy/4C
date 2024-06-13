/*! \file

\brief Factory of solid elements

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_FACTORY_HPP
#define FOUR_C_SOLID_3D_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_solid_3D_ele_calc_displacement_based.hpp"
#include "4C_solid_3D_ele_calc_displacement_based_linear_kinematics.hpp"
#include "4C_solid_3D_ele_calc_eas.hpp"
#include "4C_solid_3D_ele_calc_fbar.hpp"
#include "4C_solid_3D_ele_calc_mulf.hpp"
#include "4C_solid_3D_ele_calc_mulf_fbar.hpp"
#include "4C_solid_3D_ele_factory_lib.hpp"
#include "4C_solid_3D_ele_properties.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR::ELEMENTS
{
  enum class EasType;
}
namespace Discret::ELEMENTS
{
  namespace Details
  {
    using ImplementedSolidCellTypes = Core::FE::CelltypeSequence<Core::FE::CellType::hex8,
        Core::FE::CellType::hex18, Core::FE::CellType::hex20, Core::FE::CellType::hex27,
        Core::FE::CellType::nurbs27, Core::FE::CellType::tet4, Core::FE::CellType::tet10,
        Core::FE::CellType::wedge6, Core::FE::CellType::pyramid5>;

    using DisplacementBasedEvaluators =
        Core::FE::apply_celltype_sequence<DisplacementBasedSolidIntegrator,
            ImplementedSolidCellTypes>;

    using DisplacementBasedLinearKinematicsEvaluators =
        Core::FE::apply_celltype_sequence<DisplacementBasedLinearKinematicsSolidIntegrator,
            ImplementedSolidCellTypes>;

    using FbarEvaluators = Core::FE::apply_celltype_sequence<FBarSolidIntegrator,
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8, Core::FE::CellType::pyramid5>>;
    using EASEvaluators = Core::FE::BaseTypeList<
        SolidEleCalcEas<Core::FE::CellType::hex8, STR::ELEMENTS::EasType::eastype_h8_9>,
        SolidEleCalcEas<Core::FE::CellType::hex8, STR::ELEMENTS::EasType::eastype_h8_21>>;
    using MulfEvaluators =
        Core::FE::apply_celltype_sequence<MulfSolidIntegrator, ImplementedSolidCellTypes>;
    using FBarMulfEvaluators = Core::FE::apply_celltype_sequence<MulfFBarSolidIntegrator,
        Core::FE::CelltypeSequence<Core::FE::CellType::hex8, Core::FE::CellType::pyramid5>>;

    using SolidEvaluators =
        Core::FE::Join<DisplacementBasedEvaluators, DisplacementBasedLinearKinematicsEvaluators,
            FbarEvaluators, EASEvaluators, MulfEvaluators, FBarMulfEvaluators>;
  }  // namespace Details

  using SolidCalcVariant = CreateVariantType<Details::SolidEvaluators>;

  // forward declaration
  class SolidEleCalcInterface;
  class Solid;

  SolidCalcVariant create_solid_calculation_interface(Core::FE::CellType celltype,
      const Discret::ELEMENTS::SolidElementProperties& element_properties);

}  // namespace Discret::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif
