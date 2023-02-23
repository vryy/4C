/*----------------------------------------------------------------------*/
/*! \file
 \brief factory class providing the implementations of the porofluidmultiphase
        element evaluation routines

   \level 3

 *----------------------------------------------------------------------*/

#ifndef POROFLUIDMULTIPHASE_ELE_FACTORY_H
#define POROFLUIDMULTIPHASE_ELE_FACTORY_H


#include "lib_element.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class PoroFluidMultiPhaseEleInterface;

    class PoroFluidMultiPhaseFactory
    {
     public:
      //! ctor
      PoroFluidMultiPhaseFactory() { return; }

      //! dtor
      virtual ~PoroFluidMultiPhaseFactory() { return; }

      //! ProvideImpl
      static PoroFluidMultiPhaseEleInterface* ProvideImpl(DRT::Element::DiscretizationType distype,
          const int numdofpernode, const std::string& disname);

     private:
      //! define PoroFluidMultiPhaseEle instances dependent on problem
      template <DRT::Element::DiscretizationType distype>
      static PoroFluidMultiPhaseEleInterface* DefineProblemType(
          const int numdofpernode, const std::string& disname);


    };  // end class PoroFluidMultiPhaseFactory

  }  // namespace ELEMENTS

}  // namespace DRT



#endif  // POROFLUIDMULTIPHASE_ELE_FACTORY_H
