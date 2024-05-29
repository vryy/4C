/*! \file

\brief Factory of Shell elements

\level 1
*/

#ifndef FOUR_C_SHELL7P_ELE_FACTORY_HPP
#define FOUR_C_SHELL7P_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  // forward declaration
  class Shell7pEleCalcInterface;
  class Shell7p;

  class Shell7pFactory
  {
   public:
    //! copy constructor
    Shell7pFactory() = default;

    static std::unique_ptr<Shell7pEleCalcInterface> provide_shell7p_calculation_interface(
        const CORE::Elements::Element& ele, const std::set<INPAR::STR::EleTech>& eletech);

   private:
    //! define shell calculation instances dependent on element technology
    template <CORE::FE::CellType distype>
    static std::unique_ptr<Shell7pEleCalcInterface> define_calculation_interface_type(
        const std::set<INPAR::STR::EleTech>& eletech);
  };  // class Shell7pFactory
}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
