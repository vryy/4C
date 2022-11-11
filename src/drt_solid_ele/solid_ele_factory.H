/*----------------------------------------------------------------------*/
/*! \file

\brief Factory of solid elements

\level 1
*/
/*--------------------------------------------------------------------------*/

#ifndef SOLID_ELE_FACTORY_H
#define SOLID_ELE_FACTORY_H

#include "drt_element.H"
#include "inpar_structure.H"

namespace DRT::ELEMENTS
{
  // forward declaration
  class SolidEleInterface;
  class Solid;

  class SolidFactory
  {
   public:
    //! ctor
    SolidFactory() = default;


    static std::unique_ptr<SolidEleInterface> ProvideImpl(DRT::ELEMENTS::Solid* ele);

    template <DRT::Element::DiscretizationType distype>
    static std::unique_ptr<SolidEleInterface> ProvideImpl(DRT::ELEMENTS::Solid* ele);

  };  // class SolidFactory
}  // namespace DRT::ELEMENTS


#endif  // SOLID_ELE_FACTORY_H
