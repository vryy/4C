/*! \file
\brief Interface for scalar functions
\level 0
*/

#ifndef SRC_DRT_LIB_FUNCTION_OF_SCALAR_H
#define SRC_DRT_LIB_FUNCTION_OF_SCALAR_H

#include "function_interface_temporary.h"

namespace DRT::UTILS
{
  /**
   * This interface encodes functions \f$ y = f(s) \f$ which take a single scalar \f$ s \f$ and
   * return a single scalar \f$ y \f$.
   */
  class FunctionOfScalar : public TemporaryFunctionInterface
  {
   public:
    /**
     * Evaluate the function for the given @p scalar.
     */
    [[nodiscard]] virtual double Evaluate(double scalar) const = 0;

    /**
     * Evaluate the derivative of the function for the given @p scalar.
     */
    [[nodiscard]] virtual double EvaluateDerivative(double scalar) const = 0;
  };
}  // namespace DRT::UTILS

#endif
