/*! \file

\brief Factory of solid elements

\level 1
*/

#include "baci_solid_ele_factory.H"

#include "baci_lib_element.H"
#include "baci_solid_ele_calc.H"
#include "baci_solid_ele_calc_eas.H"
#include "baci_solid_ele_calc_fbar.H"
#include "baci_solid_ele_calc_mulf.H"
#include "baci_utils_exceptions.H"

#include <memory>

BACI_NAMESPACE_OPEN

namespace
{
  template <CORE::FE::CellType celltype>
  DRT::ELEMENTS::SolidCalcVariant CreateFBarSolidCalculationInterface(
      INPAR::STR::KinemType kinem_type)
  {
    if (kinem_type != INPAR::STR::KinemType::kinem_nonlinearTotLag)
    {
      dserror("FBAR only usable for KINEM nonlinear (you are using %s).", kinem_type);
    }

    if constexpr (celltype == CORE::FE::CellType::hex8 || celltype == CORE::FE::CellType::pyramid5)
    {
      return DRT::ELEMENTS::SolidEleCalcFbar<celltype>();
    }

    dserror("FBAR is only implemented for hex8 and pyramid5 elements.");
    return {};
  }

  template <CORE::FE::CellType celltype>
  DRT::ELEMENTS::SolidCalcVariant CreateMulfSolidCalculationInterface(
      INPAR::STR::KinemType kinem_type)
  {
    if (kinem_type != INPAR::STR::KinemType::kinem_nonlinearTotLag)
    {
      dserror("MULF only usable for KINEM nonlinear (you are using %s).", kinem_type);
    }
    return DRT::ELEMENTS::SolidEleCalcMulf<celltype>();
  }
}  // namespace

DRT::ELEMENTS::SolidCalcVariant DRT::ELEMENTS::CreateSolidCalculationInterface(
    const DRT::Element& ele, const std::set<INPAR::STR::EleTech>& eletech,
    INPAR::STR::KinemType kinem_type, STR::ELEMENTS::EasType eastype)
{
  switch (ele.Shape())
  {
    case CORE::FE::CellType::hex8:
      return CreateSolidCalculationInterface<CORE::FE::CellType::hex8>(
          ele, eletech, kinem_type, eastype);
      break;
    case CORE::FE::CellType::hex27:
      return CreateSolidCalculationInterface<CORE::FE::CellType::hex27>(
          ele, eletech, kinem_type, eastype);
      break;
    case CORE::FE::CellType::hex20:
      return CreateSolidCalculationInterface<CORE::FE::CellType::hex20>(
          ele, eletech, kinem_type, eastype);
      break;
    case CORE::FE::CellType::hex18:
      return CreateSolidCalculationInterface<CORE::FE::CellType::hex18>(
          ele, eletech, kinem_type, eastype);
      break;
    case CORE::FE::CellType::nurbs27:
      return CreateSolidCalculationInterface<CORE::FE::CellType::nurbs27>(
          ele, eletech, kinem_type, eastype);
      break;
    case CORE::FE::CellType::pyramid5:
      return CreateSolidCalculationInterface<CORE::FE::CellType::pyramid5>(
          ele, eletech, kinem_type, eastype);
      break;
    case CORE::FE::CellType::wedge6:
      return CreateSolidCalculationInterface<CORE::FE::CellType::wedge6>(
          ele, eletech, kinem_type, eastype);
      break;
    case CORE::FE::CellType::tet4:
      return CreateSolidCalculationInterface<CORE::FE::CellType::tet4>(
          ele, eletech, kinem_type, eastype);
      break;
    case CORE::FE::CellType::tet10:
      return CreateSolidCalculationInterface<CORE::FE::CellType::tet10>(
          ele, eletech, kinem_type, eastype);
      break;
    default:
      dserror("unknown celltype provided");
      break;
  }
  return {};
}

template <CORE::FE::CellType celltype>
DRT::ELEMENTS::SolidCalcVariant DRT::ELEMENTS::CreateSolidCalculationInterface(
    const DRT::Element& ele, const std::set<INPAR::STR::EleTech>& eletech,
    INPAR::STR::KinemType kinem_type, STR::ELEMENTS::EasType eastype)
{
  // here we go into the different cases for element technology
  switch (eletech.size())
  {
    // no element technology
    case 0:
      return DRT::ELEMENTS::SolidEleCalc<celltype>();
      break;
    // simple: just one element technology
    case 1:
      switch (*eletech.begin())
      {
        case INPAR::STR::EleTech::eas:
          if constexpr (celltype != CORE::FE::CellType::hex8)
          {
            dserror("EAS is only implemented for hex8 elements.");
          }
          else
          {
            switch (eastype)
            {
              case STR::ELEMENTS::EasType::eastype_h8_9:
                return DRT::ELEMENTS::SolidEleCalcEas<celltype,
                    STR::ELEMENTS::EasType::eastype_h8_9>();
              case STR::ELEMENTS::EasType::eastype_h8_21:
                return DRT::ELEMENTS::SolidEleCalcEas<celltype,
                    STR::ELEMENTS::EasType::eastype_h8_21>();
              default:
                dserror("EAS type %d is not implemented %d.", (int)eastype);
            }
          }
        case INPAR::STR::EleTech::fbar:
          return CreateFBarSolidCalculationInterface<celltype>(kinem_type);
        case INPAR::STR::EleTech::ps_mulf:
          return CreateMulfSolidCalculationInterface<celltype>(kinem_type);
        default:
          dserror("unknown element technology");
      }
    // combination of element technologies
    default:
    {
      dserror("unknown combination of element technologies.");
    }
  }
  return {};
}

BACI_NAMESPACE_CLOSE
