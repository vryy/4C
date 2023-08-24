/*! \file

\brief Factory of solid-poro elements

\level 1
*/

#ifndef BACI_SOLID_ELE_PORO_FACTORY_H
#define BACI_SOLID_ELE_PORO_FACTORY_H

#include "baci_inpar_poro.H"
#include "baci_inpar_structure.H"
#include "baci_lib_element.H"

#include <memory>

namespace DRT::ELEMENTS
{
  // forward declaration
  class SolidPoroEleCalcInterface;
  class SolidPoro;


  std::unique_ptr<SolidPoroEleCalcInterface> CreateSolidPoroCalculationInterface(
      DRT::Element& ele, INPAR::PORO::PoroType porotype);

  template <DRT::Element::DiscretizationType distype>
  std::unique_ptr<SolidPoroEleCalcInterface> CreateSolidPoroCalculationInterface(
      INPAR::PORO::PoroType porotype);
}  // namespace DRT::ELEMENTS


#endif  // SOLID_ELE_PORO_FACTORY_H