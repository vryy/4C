/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of electromagnetic elements

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_ELEMAG_ELE_FACTORY_HPP
#define FOUR_C_ELEMAG_ELE_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_lib_element.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class ElemagEleInterface;

    class ElemagFactory
    {
     public:
      //! ctor
      ElemagFactory() { return; }

      //! dtor
      virtual ~ElemagFactory() = default;
      //! ProvideImpl
      static ElemagEleInterface* ProvideImpl(CORE::FE::CellType distype, std::string problem);

     private:
      //! define ElemagEle instances dependent on problem
      template <CORE::FE::CellType distype>
      static ElemagEleInterface* DefineProblemType(std::string problem);
    };

  }  // namespace ELEMENTS

}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // ELEMAG_ELE_FACTORY_H
