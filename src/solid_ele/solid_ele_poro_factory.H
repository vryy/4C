/*----------------------------------------------------------------------*/
/*! \file

\brief Factory of solid-poro elements

\level 1
    */
/*--------------------------------------------------------------------------*/

#ifndef SOLID_ELE_PORO_FACTORY_H
#define SOLID_ELE_PORO_FACTORY_H

#include "lib_element.H"
#include "inpar_structure.H"
#include "inpar_poro.H"

namespace DRT::ELEMENTS
{
  // forward declaration
  class SolidPoroEleCalcInterface;
  class SolidPoroFactory
  {
   public:
    //! ctor
    SolidPoroFactory() = default;


    static std::unique_ptr<SolidPoroEleCalcInterface> ProvideImpl(
        DRT::Element* ele, INPAR::PORO::PoroType porotype);

    template <DRT::Element::DiscretizationType distype>
    static std::unique_ptr<SolidPoroEleCalcInterface> ProvideImpl(INPAR::PORO::PoroType porotype);

  };  // class SolidPoroFactory
}  // namespace DRT::ELEMENTS


#endif  // SOLID_ELE_PORO_FACTORY_H