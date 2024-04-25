/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of Lubrication elements

\level 3


*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_LUBRICATION_ELE_FACTORY_HPP
#define FOUR_C_LUBRICATION_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_lib_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class LubricationEleInterface;

    class LubricationFactory
    {
     public:
      //! ctor
      LubricationFactory() { return; }

      //! dtor
      virtual ~LubricationFactory() = default;
      //! ProvideImpl
      static LubricationEleInterface* ProvideImpl(
          CORE::FE::CellType distype, const std::string& disname);

     private:
      //! define LubricationEle instances dependent on problem
      template <CORE::FE::CellType distype, int probdim>
      static LubricationEleInterface* DefineProblemType(const std::string& disname);


    };  // end class LubricationFactory

  }  // namespace ELEMENTS

}  // namespace DRT


FOUR_C_NAMESPACE_CLOSE

#endif
