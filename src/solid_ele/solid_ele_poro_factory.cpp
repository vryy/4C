/*----------------------------------------------------------------------*/
/*! \file

\brief Factory of solid-poro elements

\level 1
    */
/*--------------------------------------------------------------------------*/

#include "solid_ele_poro_factory.H"
#include "solid_ele_poro_calc_pressure_based.H"
#include <memory>
#include "inpar_structure.H"



std::unique_ptr<DRT::ELEMENTS::SolidPoroEleCalcInterface>
DRT::ELEMENTS::SolidPoroFactory::ProvideImpl(DRT::Element* ele, INPAR::PORO::PoroType porotype)
{
  switch (ele->Shape())
  {
    case DRT::Element::hex8:
      return ProvideImpl<DRT::Element::hex8>(porotype);
      break;
    case DRT::Element::hex27:
      return ProvideImpl<DRT::Element::hex27>(porotype);
      break;
    case DRT::Element::hex20:
      return ProvideImpl<DRT::Element::hex20>(porotype);
      break;
    case DRT::Element::hex18:
      return ProvideImpl<DRT::Element::hex18>(porotype);
      break;
    case DRT::Element::pyramid5:
      return ProvideImpl<DRT::Element::pyramid5>(porotype);
      break;
    case DRT::Element::tet4:
      return ProvideImpl<DRT::Element::tet4>(porotype);
      break;
    case DRT::Element::tet10:
      return ProvideImpl<DRT::Element::tet10>(porotype);
      break;
    default:
      dserror("unknown distype provided");
      break;
  }
  return nullptr;
}

template <DRT::Element::DiscretizationType distype>
std::unique_ptr<DRT::ELEMENTS::SolidPoroEleCalcInterface>
DRT::ELEMENTS::SolidPoroFactory::ProvideImpl(INPAR::PORO::PoroType porotype)
{
  // here we go into the different cases for element type
  switch (porotype)
  {
    case INPAR::PORO::PoroType::porotype_pressure_velocity_based:
      dserror("POROTYPE: 'pressure_velocity_based' not yet implemented!");
      return nullptr;
      break;
    case INPAR::PORO::PoroType::porotype_pressure_based:
    {
      return std::make_unique<DRT::ELEMENTS::SolidPoroPressureBasedEleCalc<distype>>();
      break;
    }
    default:
      dserror("Wrong POROTYPE for evaluation in SolidPoro elements!");
  }
  return nullptr;
}
