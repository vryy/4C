/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of Lubrication elements

\level 3


*/
/*--------------------------------------------------------------------------*/

#ifndef LUBRICATION_ELE_FACTORY_H
#define LUBRICATION_ELE_FACTORY_H

#include "lib_element.H"

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
      virtual ~LubricationFactory() { return; }

      //! ProvideImpl
      static LubricationEleInterface* ProvideImpl(
          DRT::Element::DiscretizationType distype, const std::string& disname);

     private:
      //! define LubricationEle instances dependent on problem
      template <DRT::Element::DiscretizationType distype, int probdim>
      static LubricationEleInterface* DefineProblemType(const std::string& disname);


    };  // end class LubricationFactory

  }  // namespace ELEMENTS

}  // namespace DRT


#endif  // LUBRICATION_ELE_FACTORY_H
