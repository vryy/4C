// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_ELE_BOUNDARY_FACTORY_HPP
#define FOUR_C_FLUID_ELE_BOUNDARY_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
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

  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
