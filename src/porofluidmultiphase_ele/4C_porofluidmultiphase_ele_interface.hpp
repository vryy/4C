// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_INTERFACE_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_INTERFACE_HPP



#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace Elements
  {
    /// Interface base class for PoroFluidMultiPhaseEleCalc
    /*!
      This class exists to provide a common interface for all template
      versions of PoroFluidMultiPhaseEleCalc.
     */
    class PoroFluidMultiPhaseEleInterface
    {
     public:
      /// Empty constructor
      PoroFluidMultiPhaseEleInterface() { return; };

      /// Empty destructor
      virtual ~PoroFluidMultiPhaseEleInterface() = default;

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in PoroFluidMultiPhaseEleCalc.
       */
      virtual int evaluate(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
          std::vector<Core::LinAlg::SerialDenseMatrix*>& elemat,
          std::vector<Core::LinAlg::SerialDenseVector*>& elevec) = 0;
    };
  }  // namespace Elements
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
