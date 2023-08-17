/*----------------------------------------------------------------------*/
/*! \file

\brief Factory of solid elements

\level 1
*/
/*--------------------------------------------------------------------------*/

#include "solid_ele_factory.H"
#include <memory>
#include "utils_exceptions.H"
#include "lib_element.H"
#include "utils_exceptions.H"
#include "solid_ele_calc_eas.H"
#include "solid_ele_calc_fbar.H"
#include "solid_ele_calc.H"

namespace
{
  template <DRT::Element::DiscretizationType distype>
  std::unique_ptr<DRT::ELEMENTS::SolidEleCalcInterface> CreateFBarSolidCalculationInterface(
      INPAR::STR::KinemType kinem_type)
  {
    if constexpr (distype != DRT::Element::hex8)
    {
      dserror("FBAR is only implemented for hex8 elements.");
      return nullptr;
    }

    if (kinem_type != INPAR::STR::KinemType::kinem_nonlinearTotLag)
    {
      dserror("FBAR only usable for KINEM nonlinear (you are using %s).", kinem_type);
    }
    return std::make_unique<DRT::ELEMENTS::SolidEleCalcFbar<distype>>();
  }
}  // namespace

std::unique_ptr<DRT::ELEMENTS::SolidEleCalcInterface>
DRT::ELEMENTS::CreateSolidCalculationInterface(const DRT::Element& ele,
    const std::set<INPAR::STR::EleTech>& eletech, INPAR::STR::KinemType kinem_type,
    STR::ELEMENTS::EasType eastype)
{
  switch (ele.Shape())
  {
    case DRT::Element::hex8:
      return CreateSolidCalculationInterface<DRT::Element::hex8>(ele, eletech, kinem_type, eastype);
      break;
    case DRT::Element::hex27:
      return CreateSolidCalculationInterface<DRT::Element::hex27>(
          ele, eletech, kinem_type, eastype);
      break;
    case DRT::Element::hex20:
      return CreateSolidCalculationInterface<DRT::Element::hex20>(
          ele, eletech, kinem_type, eastype);
      break;
    case DRT::Element::hex18:
      return CreateSolidCalculationInterface<DRT::Element::hex18>(
          ele, eletech, kinem_type, eastype);
      break;
    case DRT::Element::pyramid5:
      return CreateSolidCalculationInterface<DRT::Element::pyramid5>(
          ele, eletech, kinem_type, eastype);
      break;
    case DRT::Element::wedge6:
      return CreateSolidCalculationInterface<DRT::Element::wedge6>(
          ele, eletech, kinem_type, eastype);
      break;
    case DRT::Element::tet4:
      return CreateSolidCalculationInterface<DRT::Element::tet4>(ele, eletech, kinem_type, eastype);
      break;
    case DRT::Element::tet10:
      return CreateSolidCalculationInterface<DRT::Element::tet10>(
          ele, eletech, kinem_type, eastype);
      break;
    default:
      dserror("unknown distype provided");
      break;
  }
  return nullptr;
}

template <DRT::Element::DiscretizationType distype>
std::unique_ptr<DRT::ELEMENTS::SolidEleCalcInterface>
DRT::ELEMENTS::CreateSolidCalculationInterface(const DRT::Element& ele,
    const std::set<INPAR::STR::EleTech>& eletech, INPAR::STR::KinemType kinem_type,
    STR::ELEMENTS::EasType eastype)
{
  // here we go into the different cases for element technology
  switch (eletech.size())
  {
    // no element technology
    case 0:
      return std::make_unique<DRT::ELEMENTS::SolidEleCalc<distype>>();
      break;
    // simple: just one element technology
    case 1:
      switch (*eletech.begin())
      {
        case INPAR::STR::EleTech::eas:
          if constexpr (distype != DRT::Element::hex8)
          {
            dserror("EAS is only implemented for hex8 elements.");
          }
          else
          {
            switch (eastype)
            {
              case ::STR::ELEMENTS::EasType::eastype_h8_9:
                return std::make_unique<DRT::ELEMENTS::SolidEleCalcEas<distype,
                    ::STR::ELEMENTS::EasType::eastype_h8_9>>();
              case ::STR::ELEMENTS::EasType::eastype_h8_21:
                return std::make_unique<DRT::ELEMENTS::SolidEleCalcEas<distype,
                    ::STR::ELEMENTS::EasType::eastype_h8_21>>();
              default:
                dserror("EAS type %d is not implemented %d.", (int)eastype);
            }
          }
        case INPAR::STR::EleTech::fbar:
          return CreateFBarSolidCalculationInterface<distype>(kinem_type);
        default:
          dserror("unknown element technology");
      }
    // combination of element technologies
    default:
    {
      dserror("unknown combination of element technologies.");
    }
  }
  return nullptr;
}
