/*----------------------------------------------------------------------*/
/*! \file

\brief Factory of solid elements

\level 1
*/
/*--------------------------------------------------------------------------*/

#ifndef SOLID_ELE_FACTORY_H
#define SOLID_ELE_FACTORY_H

#include "inpar_structure.H"
#include "lib_element.H"

namespace STR::ELEMENTS
{
  enum class EasType;
}
namespace DRT::ELEMENTS
{
  // forward declaration
  class SolidEleCalcInterface;
  class Solid;

  std::unique_ptr<SolidEleCalcInterface> CreateSolidCalculationInterface(const DRT::Element& ele,
      const std::set<INPAR::STR::EleTech>& eletech, INPAR::STR::KinemType kinem_type,
      STR::ELEMENTS::EasType eastype);

  template <DRT::Element::DiscretizationType distype>
  std::unique_ptr<SolidEleCalcInterface> CreateSolidCalculationInterface(const DRT::Element& ele,
      const std::set<INPAR::STR::EleTech>& eletech, INPAR::STR::KinemType kinem_type,
      STR::ELEMENTS::EasType eastype);

}  // namespace DRT::ELEMENTS


#endif  // SOLID_ELE_FACTORY_H
