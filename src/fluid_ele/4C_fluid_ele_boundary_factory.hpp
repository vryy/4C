/*----------------------------------------------------------------------*/
/*! \file

\brief factory class into templated evaluators for fluid boundary integration

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_BOUNDARY_FACTORY_HPP
#define FOUR_C_FLUID_ELE_BOUNDARY_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
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
      static FluidBoundaryInterface* provide_impl(Core::FE::CellType distype, std::string problem);

     private:
      //! define FluidEleBoundaryCalc instances dependent on problemtype
      template <Core::FE::CellType distype>
      static FluidBoundaryInterface* define_problem_type(std::string problem);
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
