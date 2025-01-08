// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_BOUNDARY_FACTORY_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    // forward declaration
    class PoroFluidMultiPhaseEleInterface;

    class PoroFluidMultiPhaseBoundaryFactory
    {
     public:
      //! ctor
      PoroFluidMultiPhaseBoundaryFactory() { return; };

      //! dtor
      virtual ~PoroFluidMultiPhaseBoundaryFactory() = default;

      //! ProvideImpl
      static PoroFluidMultiPhaseEleInterface* provide_impl(
          const Core::Elements::Element* ele, const int numdofpernode, const std::string& disname);

     private:
      //! return instance of element evaluation class depending on implementation type
      template <Core::FE::CellType distype>
      static PoroFluidMultiPhaseEleInterface* define_problem_type(
          const int numdofpernode, const std::string& disname);
    };  // class PoroFluidMultiPhaseBoundaryFactory
  }  // namespace Elements
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
