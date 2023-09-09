/*! \file

\brief Factory of Shell elements

\level 1
*/

#ifndef BACI_SHELL7P_ELE_FACTORY_H
#define BACI_SHELL7P_ELE_FACTORY_H

#include "baci_inpar_structure.H"
#include "baci_lib_element.H"

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

    static std::unique_ptr<Shell7pEleCalcInterface> ProvideShell7pCalculationInterface(
        const DRT::Element& ele, const std::set<INPAR::STR::EleTech>& eletech);

   private:
    //! define shell calculation instances dependent on element technology
    template <DRT::Element::DiscretizationType distype>
    static std::unique_ptr<Shell7pEleCalcInterface> DefineCalculationInterfaceType(
        const std::set<INPAR::STR::EleTech>& eletech);
  };  // class Shell7pFactory
}  // namespace DRT::ELEMENTS

#endif
