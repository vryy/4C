/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#include "baci_solid_ele_poro_factory.H"

#include "baci_inpar_structure.H"
#include "baci_solid_ele_poro_calc_pressure_based.H"

#include <memory>



std::unique_ptr<DRT::ELEMENTS::SolidPoroEleCalcInterface>
DRT::ELEMENTS::CreateSolidPoroCalculationInterface(
    DRT::Element& ele, INPAR::PORO::PoroType porotype)
{
  switch (ele.Shape())
  {
    case DRT::Element::hex8:
      return CreateSolidPoroCalculationInterface<DRT::Element::hex8>(porotype);
      break;
    case DRT::Element::hex27:
      return CreateSolidPoroCalculationInterface<DRT::Element::hex27>(porotype);
      break;
    case DRT::Element::hex20:
      return CreateSolidPoroCalculationInterface<DRT::Element::hex20>(porotype);
      break;
    case DRT::Element::hex18:
      return CreateSolidPoroCalculationInterface<DRT::Element::hex18>(porotype);
      break;
    case DRT::Element::pyramid5:
      return CreateSolidPoroCalculationInterface<DRT::Element::pyramid5>(porotype);
      break;
    case DRT::Element::tet4:
      return CreateSolidPoroCalculationInterface<DRT::Element::tet4>(porotype);
      break;
    case DRT::Element::tet10:
      return CreateSolidPoroCalculationInterface<DRT::Element::tet10>(porotype);
      break;
    default:
      dserror("unknown distype provided");
      break;
  }
  return nullptr;
}

template <DRT::Element::DiscretizationType distype>
std::unique_ptr<DRT::ELEMENTS::SolidPoroEleCalcInterface>
DRT::ELEMENTS::CreateSolidPoroCalculationInterface(INPAR::PORO::PoroType porotype)
{
  // here we go into the different cases for element type
  switch (porotype)
  {
    case INPAR::PORO::PoroType::pressure_velocity_based:
      dserror("POROTYPE: 'pressure_velocity_based' not yet implemented!");
      return nullptr;
      break;
    case INPAR::PORO::PoroType::pressure_based:
    {
      return std::make_unique<DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<distype>>();
      break;
    }
    default:
      dserror("Wrong POROTYPE for evaluation in SolidPoro elements!");
  }
  return nullptr;
}
