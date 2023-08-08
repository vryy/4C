/*----------------------------------------------------------------------*/
/*! \file

\brief Factory of solid elements

\level 1
*/
/*--------------------------------------------------------------------------*/

#ifndef SOLID_ELE_FACTORY_H
#define SOLID_ELE_FACTORY_H

#include "lib_element.H"
#include "inpar_structure.H"
#include "solid_ele.H"

namespace DRT::ELEMENTS
{
  // forward declaration
  class SolidEleCalcInterface;
  class Solid;

  class SolidFactory
  {
   public:
    //! ctor
    SolidFactory() = default;


    static std::unique_ptr<SolidEleCalcInterface> ProvideImpl(DRT::Element* ele,
        std::set<INPAR::STR::EleTech> eletech, ::INPAR::STR::KinemType kinemtype,
        ::STR::ELEMENTS::EASType eastype);

    template <DRT::Element::DiscretizationType distype>
    static std::unique_ptr<SolidEleCalcInterface> ProvideImpl(std::set<INPAR::STR::EleTech> eletech,
        ::INPAR::STR::KinemType kinemtype, ::STR::ELEMENTS::EASType eastype);

  };  // class SolidFactory
}  // namespace DRT::ELEMENTS


#endif  // SOLID_ELE_FACTORY_H
