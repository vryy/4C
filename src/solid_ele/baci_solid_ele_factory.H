/*! \file

\brief Factory of solid elements

\level 1
*/

#ifndef BACI_SOLID_ELE_FACTORY_H
#define BACI_SOLID_ELE_FACTORY_H

#include "baci_config.H"

#include "baci_discretization_fem_general_cell_type_traits.H"
#include "baci_solid_ele_calc_displacement_based.H"
#include "baci_solid_ele_calc_eas.H"
#include "baci_solid_ele_calc_fbar.H"
#include "baci_solid_ele_calc_mulf.H"
#include "baci_solid_ele_factory_lib.H"
#include "baci_solid_ele_properties.H"

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

    using FbarEvaluators = CORE::FE::apply_celltype_sequence<SolidEleCalcFbar,
        CORE::FE::celltype_sequence<CORE::FE::CellType::hex8, CORE::FE::CellType::pyramid5>>;
    using EASEvaluators = CORE::FE::BaseTypeList<
        SolidEleCalcEas<CORE::FE::CellType::hex8, STR::ELEMENTS::EasType::eastype_h8_9>,
        SolidEleCalcEas<CORE::FE::CellType::hex8, STR::ELEMENTS::EasType::eastype_h8_21>>;
    using MulfEvaluators =
        CORE::FE::apply_celltype_sequence<SolidEleCalcMulf, ImplementedSolidCellTypes>;

    using SolidEvaluators =
        CORE::FE::Join<DisplacementBasedEvaluators, FbarEvaluators, EASEvaluators, MulfEvaluators>;
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
