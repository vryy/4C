/*----------------------------------------------------------------------*/
/*! \file

\brief factory class into templated evaluators for fluid boundary integration

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_BOUNDARY_FACTORY_H
#define BACI_FLUID_ELE_BOUNDARY_FACTORY_H

#include "baci_lib_element.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidBoundaryInterface;

    class FluidBoundaryFactory
    {
     public:
      //! ctor
      FluidBoundaryFactory() { return; };

      //! dtor
      virtual ~FluidBoundaryFactory() { return; };

      //! ProvideImpl
      static FluidBoundaryInterface* ProvideImpl(
          DRT::Element::DiscretizationType distype, std::string problem);

     private:
      //! define FluidEleBoundaryCalc instances dependent on problemtype
      template <DRT::Element::DiscretizationType distype>
      static FluidBoundaryInterface* DefineProblemType(std::string problem);
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
