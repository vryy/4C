/*! \file

\brief Factory of Shell7p elements

\level 1
*/


#include "baci_shell7p_ele_factory.H"

#include "baci_shell7p_ele.H"
#include "baci_shell7p_ele_calc.H"
#include "baci_shell7p_ele_calc_eas.H"

#include <memory>

std::unique_ptr<DRT::ELEMENTS::Shell7pEleCalcInterface>
DRT::ELEMENTS::Shell7pFactory::ProvideShell7pCalculationInterface(
    const DRT::Element& ele, const std::set<INPAR::STR::EleTech>& eletech)
{
  switch (ele.Shape())
  {
    case DRT::Element::quad4:
    {
      return DefineCalculationInterfaceType<DRT::Element::quad4>(eletech);
    }
    case DRT::Element::quad8:
    {
      return DefineCalculationInterfaceType<DRT::Element::quad8>(eletech);
    }
    case DRT::Element::quad9:
    {
      return DefineCalculationInterfaceType<DRT::Element::quad9>(eletech);
    }
    case DRT::Element::tri3:
    {
      return DefineCalculationInterfaceType<DRT::Element::tri3>(eletech);
    }
    case DRT::Element::tri6:
    {
      return DefineCalculationInterfaceType<DRT::Element::tri6>(eletech);
    }
    default:
      dserror("unknown distype provided");
  }
}

template <DRT::Element::DiscretizationType distype>
std::unique_ptr<DRT::ELEMENTS::Shell7pEleCalcInterface>
DRT::ELEMENTS::Shell7pFactory::DefineCalculationInterfaceType(
    const std::set<INPAR::STR::EleTech>& eletech)
{
  // here we go into the different cases for element technology
  switch (eletech.size())
  {
    // no element technology
    case 0:
      return std::make_unique<DRT::ELEMENTS::Shell7pEleCalc<distype>>();
    // simple: just one element technology
    case 1:
      switch (*eletech.begin())
      {
        case INPAR::STR::EleTech::eas:
        {
          if constexpr ((distype != DRT::Element::quad4) && (distype != DRT::Element::quad9))
          {
            dserror("EAS is only implemented for quad4 and quad9 elements.");
          }
          return std::make_unique<DRT::ELEMENTS::Shell7pEleCalcEas<distype>>();
        }
        default:
          dserror("unknown element technology");
      }
    // combination of element technologies
    default:
    {
      dserror("unknown combination of element technologies.");
    }
  }
}
