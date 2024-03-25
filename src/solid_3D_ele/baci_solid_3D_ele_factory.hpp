/*! \file

\brief Factory of solid elements

\level 1
*/

#ifndef BACI_SOLID_3D_ELE_FACTORY_HPP
#define BACI_SOLID_3D_ELE_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_solid_3D_ele_calc_displacement_based.hpp"
#include "baci_solid_3D_ele_calc_displacement_based_linear_kinematics.hpp"
#include "baci_solid_3D_ele_calc_eas.hpp"
#include "baci_solid_3D_ele_calc_fbar.hpp"
#include "baci_solid_3D_ele_calc_mulf.hpp"
#include "baci_solid_3D_ele_calc_mulf_fbar.hpp"
#include "baci_solid_3D_ele_factory_lib.hpp"
#include "baci_solid_3D_ele_properties.hpp"

BACI_NAMESPACE_OPEN

namespace STR::ELEMENTS
{
  enum class EasType;
}
namespace DRT::ELEMENTS
{
  namespace DETAILS
  {
    using ImplementedSolidCellTypes = CORE::FE::celltype_sequence<CORE::FE::CellType::hex8,
        CORE::FE::CellType::hex18, CORE::FE::CellType::hex20, CORE::FE::CellType::hex27,
        CORE::FE::CellType::nurbs27, CORE::FE::CellType::tet4, CORE::FE::CellType::tet10,
        CORE::FE::CellType::wedge6, CORE::FE::CellType::pyramid5>;

    using DisplacementBasedEvaluators =
        CORE::FE::apply_celltype_sequence<DisplacementBasedSolidIntegrator,
            ImplementedSolidCellTypes>;

    using DisplacementBasedLinearKinematicsEvaluators =
        CORE::FE::apply_celltype_sequence<DisplacementBasedLinearKinematicsSolidIntegrator,
            ImplementedSolidCellTypes>;

    using FbarEvaluators = CORE::FE::apply_celltype_sequence<FBarSolidIntegrator,
        CORE::FE::celltype_sequence<CORE::FE::CellType::hex8, CORE::FE::CellType::pyramid5>>;
    using EASEvaluators = CORE::FE::BaseTypeList<
        SolidEleCalcEas<CORE::FE::CellType::hex8, STR::ELEMENTS::EasType::eastype_h8_9>,
        SolidEleCalcEas<CORE::FE::CellType::hex8, STR::ELEMENTS::EasType::eastype_h8_21>>;
    using MulfEvaluators =
        CORE::FE::apply_celltype_sequence<MulfSolidIntegrator, ImplementedSolidCellTypes>;
    using FBarMulfEvaluators = CORE::FE::apply_celltype_sequence<MulfFBarSolidIntegrator,
        CORE::FE::celltype_sequence<CORE::FE::CellType::hex8, CORE::FE::CellType::pyramid5>>;

    using SolidEvaluators =
        CORE::FE::Join<DisplacementBasedEvaluators, DisplacementBasedLinearKinematicsEvaluators,
            FbarEvaluators, EASEvaluators, MulfEvaluators, FBarMulfEvaluators>;
  }  // namespace DETAILS

  using SolidCalcVariant = CreateVariantType<DETAILS::SolidEvaluators>;

  // forward declaration
  class SolidEleCalcInterface;
  class Solid;

  SolidCalcVariant CreateSolidCalculationInterface(
      CORE::FE::CellType celltype, const DRT::ELEMENTS::SolidElementProperties& element_properties);

}  // namespace DRT::ELEMENTS


BACI_NAMESPACE_CLOSE

#endif  // SOLID_ELE_FACTORY_H
