/*! \file

\brief Factory of Shell7p elements

\level 1
*/


#include "4C_shell7p_ele_factory.hpp"

#include "4C_shell7p_ele.hpp"
#include "4C_shell7p_ele_calc.hpp"
#include "4C_shell7p_ele_calc_eas.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

std::unique_ptr<Discret::ELEMENTS::Shell7pEleCalcInterface>
Discret::ELEMENTS::Shell7pFactory::provide_shell7p_calculation_interface(
    const Core::Elements::Element& ele, const std::set<Inpar::STR::EleTech>& eletech)
{
  switch (ele.Shape())
  {
    case Core::FE::CellType::quad4:
    {
      return define_calculation_interface_type<Core::FE::CellType::quad4>(eletech);
    }
    case Core::FE::CellType::quad8:
    {
      return define_calculation_interface_type<Core::FE::CellType::quad8>(eletech);
    }
    case Core::FE::CellType::quad9:
    {
      return define_calculation_interface_type<Core::FE::CellType::quad9>(eletech);
    }
    case Core::FE::CellType::tri3:
    {
      return define_calculation_interface_type<Core::FE::CellType::tri3>(eletech);
    }
    case Core::FE::CellType::tri6:
    {
      return define_calculation_interface_type<Core::FE::CellType::tri6>(eletech);
    }
    default:
      FOUR_C_THROW("unknown distype provided");
  }
}

template <Core::FE::CellType distype>
std::unique_ptr<Discret::ELEMENTS::Shell7pEleCalcInterface>
Discret::ELEMENTS::Shell7pFactory::define_calculation_interface_type(
    const std::set<Inpar::STR::EleTech>& eletech)
{
  // here we go into the different cases for element technology
  switch (eletech.size())
  {
    // no element technology
    case 0:
      return std::make_unique<Discret::ELEMENTS::Shell7pEleCalc<distype>>();
    // simple: just one element technology
    case 1:
      switch (*eletech.begin())
      {
        case Inpar::STR::EleTech::eas:
        {
          if constexpr ((distype != Core::FE::CellType::quad4) &&
                        (distype != Core::FE::CellType::quad9))
          {
            FOUR_C_THROW("EAS is only implemented for quad4 and quad9 elements.");
          }
          return std::make_unique<Discret::ELEMENTS::Shell7pEleCalcEas<distype>>();
        }
        default:
          FOUR_C_THROW("unknown element technology");
      }
    // combination of element technologies
    default:
    {
      FOUR_C_THROW("unknown combination of element technologies.");
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
