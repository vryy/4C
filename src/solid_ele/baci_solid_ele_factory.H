/*! \file

\brief Factory of solid elements

\level 1
*/

#ifndef BACI_SOLID_ELE_FACTORY_H
#define BACI_SOLID_ELE_FACTORY_H

#include "baci_discretization_fem_general_cell_type_traits.H"
#include "baci_inpar_structure.H"
#include "baci_lib_element.H"
#include "baci_solid_ele_calc.H"
#include "baci_solid_ele_calc_eas.H"
#include "baci_solid_ele_calc_fbar.H"
#include "baci_solid_ele_calc_mulf.H"
#include "baci_solid_ele_factory_lib.H"

#include <variant>

namespace STR::ELEMENTS
{
  enum class EasType;
}
namespace DRT::ELEMENTS
{
  /*!
   *  @brief struct for managing solid element properties
   */
  struct SolidElementProperties
  {
    //! kinematic type
    INPAR::STR::KinemType kintype{INPAR::STR::kinem_vague};

    //! element technology
    std::set<INPAR::STR::EleTech> eletech;

    //! specify EAS type in case contains enhanced assumed strains
    STR::ELEMENTS::EasType eastype{STR::ELEMENTS::EasType::eastype_undefined};
  };

  namespace DETAILS
  {
    using ImplementedSolidCellTypes = CORE::FE::celltype_sequence<CORE::FE::CellType::hex8,
        CORE::FE::CellType::hex18, CORE::FE::CellType::hex20, CORE::FE::CellType::hex27,
        CORE::FE::CellType::nurbs27, CORE::FE::CellType::tet4, CORE::FE::CellType::tet10,
        CORE::FE::CellType::wedge6, CORE::FE::CellType::pyramid5>;

    using DisplacementBasedEvaluators =
        CORE::FE::apply_celltype_sequence<SolidEleCalc, ImplementedSolidCellTypes>;

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

  SolidCalcVariant CreateSolidCalculationInterface(const DRT::Element& ele,
      const std::set<INPAR::STR::EleTech>& eletech, INPAR::STR::KinemType kinem_type,
      STR::ELEMENTS::EasType eastype);

  template <CORE::FE::CellType distype>
  SolidCalcVariant CreateSolidCalculationInterface(const DRT::Element& ele,
      const std::set<INPAR::STR::EleTech>& eletech, INPAR::STR::KinemType kinem_type,
      STR::ELEMENTS::EasType eastype);

}  // namespace DRT::ELEMENTS


#endif  // SOLID_ELE_FACTORY_H
