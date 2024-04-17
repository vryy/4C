/*----------------------------------------------------------------------*/
/*! \file

\brief factory class into templated evaluators for fluid boundary integration

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_BOUNDARY_FACTORY_HPP
#define FOUR_C_FLUID_ELE_BOUNDARY_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_lib_element.hpp"

FOUR_C_NAMESPACE_OPEN

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
      virtual ~FluidBoundaryFactory() = default;

      //! ProvideImpl
      static FluidBoundaryInterface* ProvideImpl(CORE::FE::CellType distype, std::string problem);

     private:
      //! define FluidEleBoundaryCalc instances dependent on problemtype
      template <CORE::FE::CellType distype>
      static FluidBoundaryInterface* DefineProblemType(std::string problem);
    };

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
