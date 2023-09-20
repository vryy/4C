/*! \file

\brief Factory of solid elements

\level 1
*/

#ifndef BACI_SOLID_ELE_FACTORY_H
#define BACI_SOLID_ELE_FACTORY_H

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
  namespace DETAILS
  {
    using ImplementedSolidDiscretizationTypes =
        DRT::DiscretizationTypeList<DRT::Element::DiscretizationType::hex8,
            DRT::Element::DiscretizationType::hex18, DRT::Element::DiscretizationType::hex20,
            DRT::Element::DiscretizationType::hex27, DRT::Element::DiscretizationType::tet4,
            DRT::Element::DiscretizationType::tet10, DRT::Element::DiscretizationType::wedge6,
            DRT::Element::DiscretizationType::pyramid5>;

    using DisplacementBasedEvaluators =
        DRT::ApplyToBase<SolidEleCalc, ImplementedSolidDiscretizationTypes>;
    using FbarEvaluators = DRT::ApplyToBase<SolidEleCalcFbar,
        DRT::DiscretizationTypeList<DRT::Element::DiscretizationType::hex8,
            DRT::Element::DiscretizationType::pyramid5>>;
    using EASEvaluators = DRT::BaseTypeList<SolidEleCalcEas<DRT::Element::DiscretizationType::hex8,
                                                STR::ELEMENTS::EasType::eastype_h8_9>,
        SolidEleCalcEas<DRT::Element::DiscretizationType::hex8,
            STR::ELEMENTS::EasType::eastype_h8_21>>;
    using MulfEvaluators = DRT::ApplyToBase<SolidEleCalcMulf, ImplementedSolidDiscretizationTypes>;

    using SolidEvaluators =
        DRT::Join_t<DisplacementBasedEvaluators, FbarEvaluators, EASEvaluators, MulfEvaluators>;
  }  // namespace DETAILS

  using SolidCalcVariant = CreateVariantType<DETAILS::SolidEvaluators>;

  // forward declaration
  class SolidEleCalcInterface;
  class Solid;

  SolidCalcVariant CreateSolidCalculationInterface(const DRT::Element& ele,
      const std::set<INPAR::STR::EleTech>& eletech, INPAR::STR::KinemType kinem_type,
      STR::ELEMENTS::EasType eastype);

  template <DRT::Element::DiscretizationType distype>
  SolidCalcVariant CreateSolidCalculationInterface(const DRT::Element& ele,
      const std::set<INPAR::STR::EleTech>& eletech, INPAR::STR::KinemType kinem_type,
      STR::ELEMENTS::EasType eastype);

}  // namespace DRT::ELEMENTS


#endif  // SOLID_ELE_FACTORY_H
