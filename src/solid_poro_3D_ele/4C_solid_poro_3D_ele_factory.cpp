/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#include "4C_solid_poro_3D_ele_factory.hpp"

#include "4C_solid_poro_3D_ele_calc_pressure_based.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN



Discret::ELEMENTS::SolidPoroCalcVariant Discret::ELEMENTS::CreateSolidPoroCalculationInterface(
    Core::Elements::Element& ele, Inpar::Poro::PoroType porotype)
{
  switch (ele.Shape())
  {
    case Core::FE::CellType::hex8:
      return CreateSolidPoroCalculationInterface<Core::FE::CellType::hex8>(porotype);
      break;
    case Core::FE::CellType::hex27:
      return CreateSolidPoroCalculationInterface<Core::FE::CellType::hex27>(porotype);
      break;
    case Core::FE::CellType::tet4:
      return CreateSolidPoroCalculationInterface<Core::FE::CellType::tet4>(porotype);
      break;
    case Core::FE::CellType::tet10:
      return CreateSolidPoroCalculationInterface<Core::FE::CellType::tet10>(porotype);
      break;
    default:
      FOUR_C_THROW("unknown celltype provided");
      break;
  }
  return {};
}

template <Core::FE::CellType celltype>
Discret::ELEMENTS::SolidPoroCalcVariant Discret::ELEMENTS::CreateSolidPoroCalculationInterface(
    Inpar::Poro::PoroType porotype)
{
  // here we go into the different cases for poro type
  switch (porotype)
  {
    case Inpar::Poro::PoroType::pressure_velocity_based:
      FOUR_C_THROW("POROTYPE: 'pressure_velocity_based' not yet implemented!");
      return {};
      break;
    case Inpar::Poro::PoroType::pressure_based:
    {
      return Discret::ELEMENTS::SolidPoroPressureBasedEleCalc<celltype>();
      break;
    }
    default:
      FOUR_C_THROW("Wrong POROTYPE for evaluation in SolidPoro elements!");
  }
  return {};
}

FOUR_C_NAMESPACE_CLOSE
