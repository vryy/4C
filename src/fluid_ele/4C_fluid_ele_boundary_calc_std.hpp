/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of fluid terms at integration points of boundaries

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_BOUNDARY_CALC_STD_HPP
#define FOUR_C_FLUID_ELE_BOUNDARY_CALC_STD_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_boundary_calc.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class FluidBoundary;

    template <Core::FE::CellType distype>
    class FluidEleBoundaryCalcStd : public FluidBoundaryImpl<distype>
    {
      typedef Discret::ELEMENTS::FluidBoundaryImpl<distype> my;

     public:
      /// Singleton access method
      static FluidEleBoundaryCalcStd<distype>* Instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

     private:
      /// private constructor since we are singleton
      FluidEleBoundaryCalcStd();

    };  // class FluidEleBoundaryCalcStd

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
