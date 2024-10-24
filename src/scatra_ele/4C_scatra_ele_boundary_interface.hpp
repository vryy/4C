// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_INTERFACE_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class FaceElement;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace Elements
  {
    class TransportBoundary;

    /// Interface base class for ScaTraEleBoundaryCalc
    /*!
      This class exists to provide a common interface for all template
      versions of ScaTraEleBoundaryCalc.
     */
    class ScaTraBoundaryInterface
    {
     public:
      /// Virtual destructor.
      virtual ~ScaTraBoundaryInterface() = default;

      virtual int evaluate(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) = 0;

      virtual int evaluate_neumann(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
          Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseVector& elevec1,
          const double scalar) = 0;
    };
  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
