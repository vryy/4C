/*! \file

\brief Factory of Shell7p elements

\level 1
*/


#include "baci_shell7p_ele_factory.hpp"

#include "baci_shell7p_ele.hpp"
#include "baci_shell7p_ele_calc.hpp"
#include "baci_shell7p_ele_calc_eas.hpp"

#include <memory>

BACI_NAMESPACE_OPEN

std::unique_ptr<DRT::ELEMENTS::Shell7pEleCalcInterface>
DRT::ELEMENTS::Shell7pFactory::ProvideShell7pCalculationInterface(
    const DRT::Element& ele, const std::set<INPAR::STR::EleTech>& eletech)
{
  switch (ele.Shape())
  {
    case CORE::FE::CellType::quad4:
    {
      return DefineCalculationInterfaceType<CORE::FE::CellType::quad4>(eletech);
    }
    case CORE::FE::CellType::quad8:
    {
      return DefineCalculationInterfaceType<CORE::FE::CellType::quad8>(eletech);
    }
    case CORE::FE::CellType::quad9:
    {
      return DefineCalculationInterfaceType<CORE::FE::CellType::quad9>(eletech);
    }
    case CORE::FE::CellType::tri3:
    {
      return DefineCalculationInterfaceType<CORE::FE::CellType::tri3>(eletech);
    }
    case CORE::FE::CellType::tri6:
    {
      return DefineCalculationInterfaceType<CORE::FE::CellType::tri6>(eletech);
    }
    default:
      dserror("unknown distype provided");
  }
}

template <CORE::FE::CellType distype>
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
          if constexpr ((distype != CORE::FE::CellType::quad4) &&
                        (distype != CORE::FE::CellType::quad9))
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

BACI_NAMESPACE_CLOSE
