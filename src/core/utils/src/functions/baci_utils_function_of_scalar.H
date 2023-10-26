/*! \file
\brief Interface for scalar functions
\level 0
*/

#ifndef BACI_UTILS_FUNCTION_OF_SCALAR_H
#define BACI_UTILS_FUNCTION_OF_SCALAR_H

namespace CORE::UTILS
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
    [[nodiscard]] virtual double Evaluate(double scalar) const = 0;

    /**
     * Evaluate the @deriv_order derivative of the function for the given @p scalar.
     */
    [[nodiscard]] virtual double EvaluateDerivative(double scalar, int deriv_order) const = 0;
  };
}  // namespace CORE::UTILS

#endif
