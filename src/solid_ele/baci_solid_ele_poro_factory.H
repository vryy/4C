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
        DRT::DiscretizationTypeList<DRT::Element::DiscretizationType::hex8,
            DRT::Element::DiscretizationType::hex18, DRT::Element::DiscretizationType::hex20,
            DRT::Element::DiscretizationType::hex27, DRT::Element::DiscretizationType::tet4,
            DRT::Element::DiscretizationType::tet10, DRT::Element::DiscretizationType::pyramid5>;
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

  template <DRT::Element::DiscretizationType distype>
  SolidPoroCalcVariant CreateSolidPoroCalculationInterface(INPAR::PORO::PoroType porotype);
}  // namespace DRT::ELEMENTS


#endif  // SOLID_ELE_PORO_FACTORY_H