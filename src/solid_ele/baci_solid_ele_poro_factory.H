/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#ifndef BACI_SOLID_ELE_PORO_FACTORY_H
#define BACI_SOLID_ELE_PORO_FACTORY_H

#include "baci_inpar_poro.H"
#include "baci_lib_element.H"
#include "baci_solid_ele_factory_lib.H"
#include "baci_solid_ele_poro_calc_pressure_based.H"

#include <memory>
#include <variant>

namespace DRT::ELEMENTS
{
  namespace DETAILS
  {
    using ImplementedSolidPoroDiscretizationTypes =
        DRT::DiscretizationTypeList<CORE::FE::CellType::hex8, CORE::FE::CellType::hex18,
            CORE::FE::CellType::hex20, CORE::FE::CellType::hex27, CORE::FE::CellType::tet4,
            CORE::FE::CellType::tet10, CORE::FE::CellType::pyramid5>;
    using PoroPressureBasedEvaluators =
        DRT::ApplyToBase<SolidPoroPressureBasedEleCalc, ImplementedSolidPoroDiscretizationTypes>;

    using SolidPoroEvaluators = DRT::Join_t<PoroPressureBasedEvaluators>;
  }  // namespace DETAILS
  using SolidPoroCalcVariant = CreateVariantType<DETAILS::SolidPoroEvaluators>;

  // forward declaration
  class SolidPoroEleCalcInterface;
  class SolidPoro;


  SolidPoroCalcVariant CreateSolidPoroCalculationInterface(
      DRT::Element& ele, INPAR::PORO::PoroType porotype);

  template <CORE::FE::CellType distype>
  SolidPoroCalcVariant CreateSolidPoroCalculationInterface(INPAR::PORO::PoroType porotype);
}  // namespace DRT::ELEMENTS


#endif  // SOLID_ELE_PORO_FACTORY_H