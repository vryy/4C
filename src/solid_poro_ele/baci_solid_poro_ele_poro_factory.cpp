/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#include "baci_solid_poro_ele_calc_pressure_based.H"
#include "baci_solid_poro_ele_factory.H"

#include <memory>



DRT::ELEMENTS::SolidPoroCalcVariant DRT::ELEMENTS::CreateSolidPoroCalculationInterface(
    DRT::Element& ele, INPAR::PORO::PoroType porotype)
{
  switch (ele.Shape())
  {
    case CORE::FE::CellType::hex8:
      return CreateSolidPoroCalculationInterface<CORE::FE::CellType::hex8>(porotype);
      break;
    case CORE::FE::CellType::hex27:
      return CreateSolidPoroCalculationInterface<CORE::FE::CellType::hex27>(porotype);
      break;
    case CORE::FE::CellType::hex20:
      return CreateSolidPoroCalculationInterface<CORE::FE::CellType::hex20>(porotype);
      break;
    case CORE::FE::CellType::hex18:
      return CreateSolidPoroCalculationInterface<CORE::FE::CellType::hex18>(porotype);
      break;
    case CORE::FE::CellType::pyramid5:
      return CreateSolidPoroCalculationInterface<CORE::FE::CellType::pyramid5>(porotype);
      break;
    case CORE::FE::CellType::tet4:
      return CreateSolidPoroCalculationInterface<CORE::FE::CellType::tet4>(porotype);
      break;
    case CORE::FE::CellType::tet10:
      return CreateSolidPoroCalculationInterface<CORE::FE::CellType::tet10>(porotype);
      break;
    default:
      dserror("unknown distype provided");
      break;
  }
  return {};
}

template <CORE::FE::CellType distype>
DRT::ELEMENTS::SolidPoroCalcVariant DRT::ELEMENTS::CreateSolidPoroCalculationInterface(
    INPAR::PORO::PoroType porotype)
{
  // here we go into the different cases for poro type
  switch (porotype)
  {
    case INPAR::PORO::PoroType::pressure_velocity_based:
      dserror("POROTYPE: 'pressure_velocity_based' not yet implemented!");
      return {};
      break;
    case INPAR::PORO::PoroType::pressure_based:
    {
      return DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<distype>();
      break;
    }
    default:
      dserror("Wrong POROTYPE for evaluation in SolidPoro elements!");
  }
  return {};
}
