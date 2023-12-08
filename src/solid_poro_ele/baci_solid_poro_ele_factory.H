/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#ifndef BACI_SOLID_PORO_ELE_FACTORY_H
#define BACI_SOLID_PORO_ELE_FACTORY_H

#include "baci_discretization_fem_general_cell_type_traits.H"
#include "baci_inpar_poro.H"
#include "baci_inpar_scatra.H"
#include "baci_lib_element.H"
#include "baci_solid_ele_factory_lib.H"
#include "baci_solid_poro_ele_calc_pressure_based.H"

#include <memory>
#include <variant>

namespace DRT::ELEMENTS
{
  /*!
   *  @brief struct for managing solidporo element properties
   */
  struct SolidPoroElementProperties
  {
    //! porosity implementation type (physics)
    INPAR::PORO::PoroType porotype{INPAR::PORO::PoroType::undefined};

    //! scalar transport implementation type (physics)
    INPAR::SCATRA::ImplType impltype{INPAR::SCATRA::ImplType::impltype_undefined};
  };

  namespace DETAILS
  {
    using ImplementedSolidPoroCellTypes = CORE::FE::celltype_sequence<CORE::FE::CellType::hex8,
        CORE::FE::CellType::hex27, CORE::FE::CellType::tet4, CORE::FE::CellType::tet10>;
    using PoroPressureBasedEvaluators =
        CORE::FE::apply_celltype_sequence<SolidPoroPressureBasedEleCalc,
            ImplementedSolidPoroCellTypes>;

    using SolidPoroEvaluators = CORE::FE::Join<PoroPressureBasedEvaluators>;
  }  // namespace DETAILS
  using SolidPoroCalcVariant = CreateVariantType<DETAILS::SolidPoroEvaluators>;

  // forward declaration
  class SolidPoroEleCalcInterface;
  class SolidPoro;


  SolidPoroCalcVariant CreateSolidPoroCalculationInterface(
      DRT::Element& ele, INPAR::PORO::PoroType porotype);

  template <CORE::FE::CellType celltype>
  SolidPoroCalcVariant CreateSolidPoroCalculationInterface(INPAR::PORO::PoroType porotype);
}  // namespace DRT::ELEMENTS


#endif  // BACI_SOLID_PORO_ELE_FACTORY_H