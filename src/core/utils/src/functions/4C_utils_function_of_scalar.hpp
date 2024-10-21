// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_FUNCTION_OF_SCALAR_HPP
#define FOUR_C_UTILS_FUNCTION_OF_SCALAR_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /**
   * This interface encodes functions \f$ y = f(s) \f$ which take a single scalar \f$ s \f$ and
   * return a single scalar \f$ y \f$.
   */
  class FunctionOfScalar
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~FunctionOfScalar() = default;

    /**
     * Evaluate the function for the given @p scalar.
     */
    [[nodiscard]] virtual double evaluate(double scalar) const = 0;

    /**
     * Evaluate the @deriv_order derivative of the function for the given @p scalar.
     */
    [[nodiscard]] virtual double evaluate_derivative(double scalar, int deriv_order) const = 0;
  };
}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
