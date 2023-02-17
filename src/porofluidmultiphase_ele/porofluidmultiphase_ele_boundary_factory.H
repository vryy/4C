/*----------------------------------------------------------------------*/
/*! \file
 \brief  factory class providing the implementations of the porous fluid multiphase
         boundary element

   \level 3

 *----------------------------------------------------------------------*/

#ifndef POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_H
#define POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_H


#include "lib_element.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class PoroFluidMultiPhaseEleInterface;

    class PoroFluidMultiPhaseBoundaryFactory
    {
     public:
      //! ctor
      PoroFluidMultiPhaseBoundaryFactory() { return; };

      //! dtor
      virtual ~PoroFluidMultiPhaseBoundaryFactory() { return; };

      //! ProvideImpl
      static PoroFluidMultiPhaseEleInterface* ProvideImpl(
          const DRT::Element* ele, const int numdofpernode, const std::string& disname);

     private:
      //! return instance of element evaluation class depending on implementation type
      template <DRT::Element::DiscretizationType distype>
      static PoroFluidMultiPhaseEleInterface* DefineProblemType(
          const int numdofpernode, const std::string& disname);
    };  // class PoroFluidMultiPhaseBoundaryFactory
  }     // namespace ELEMENTS
}  // namespace DRT


#endif  // POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_H
