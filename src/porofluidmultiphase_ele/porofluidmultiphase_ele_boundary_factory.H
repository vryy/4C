/*----------------------------------------------------------------------*/
/*! \file
 \brief  factory class providing the implementations of the porous fluid multiphase
         boundary element

   \level 3

 *----------------------------------------------------------------------*/

#ifndef SRC_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_H_
#define SRC_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_H_


#include "element.H"

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


#endif /* SRC_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_H_ */
