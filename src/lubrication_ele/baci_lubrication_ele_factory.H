/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of Lubrication elements

\level 3


*/
/*--------------------------------------------------------------------------*/

#ifndef BACI_LUBRICATION_ELE_FACTORY_H
#define BACI_LUBRICATION_ELE_FACTORY_H

#include "baci_lib_element.H"

BACI_NAMESPACE_OPEN

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


BACI_NAMESPACE_CLOSE

#endif  // LUBRICATION_ELE_FACTORY_H
