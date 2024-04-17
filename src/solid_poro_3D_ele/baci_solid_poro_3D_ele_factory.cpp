/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#include "baci_solid_poro_3D_ele_factory.hpp"

#include "baci_solid_poro_3D_ele_calc_pressure_based.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN



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
    case CORE::FE::CellType::tet4:
      return CreateSolidPoroCalculationInterface<CORE::FE::CellType::tet4>(porotype);
      break;
    case CORE::FE::CellType::tet10:
      return CreateSolidPoroCalculationInterface<CORE::FE::CellType::tet10>(porotype);
      break;
    default:
      dserror("unknown celltype provided");
      break;
  }
  return {};
}

template <CORE::FE::CellType celltype>
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
      return DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<celltype>();
      break;
    }
    default:
      dserror("Wrong POROTYPE for evaluation in SolidPoro elements!");
  }
  return {};
}

FOUR_C_NAMESPACE_CLOSE
