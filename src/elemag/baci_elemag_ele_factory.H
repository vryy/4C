/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of electromagnetic elements

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef BACI_ELEMAG_ELE_FACTORY_H
#define BACI_ELEMAG_ELE_FACTORY_H

#include "baci_lib_element.H"

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

#endif  // ELEMAG_ELE_FACTORY_H
