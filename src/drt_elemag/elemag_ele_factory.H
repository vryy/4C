/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of electromagnetic elements

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef ELEMAG_ELE_FACTORY_H
#define ELEMAG_ELE_FACTORY_H

#include "../drt_lib/drt_element.H"

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
      virtual ~ElemagFactory() { return; }

      //! ProvideImpl
      static ElemagEleInterface* ProvideImpl(
          DRT::Element::DiscretizationType distype, std::string problem);

     private:
      //! define ElemagEle instances dependent on problem
      template <DRT::Element::DiscretizationType distype>
      static ElemagEleInterface* DefineProblemType(std::string problem);
    };

  }  // namespace ELEMENTS

}  // namespace DRT

#endif /* ELEMAG_ELE_FACTORY_H_ */
