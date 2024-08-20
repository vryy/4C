/*! \file
\brief Interface for functions of time
\level 0
*/

#ifndef FOUR_C_UTILS_FUNCTION_OF_TIME_HPP
#define FOUR_C_UTILS_FUNCTION_OF_TIME_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_polynomial.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function_manager.hpp"
#include "4C_utils_functionvariables.hpp"
#include "4C_utils_symbolic_expression.hpp"

#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

#include <complex>
#include <iostream>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::UTILS
{
  /*!
   * \brief interface for time-dependent functions.
   *
   * It encodes potentially vector-valued functions \f$ y_i = f_i(t) \f$ which take a time value
   * \f$ t \f$ and return the component \f$ y_i \f$ or its first derivative.

   */
  class FunctionOfTime
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~FunctionOfTime() = default;

    /**
     * Evaluate the function for the given @p time and @p component.
     */
    [[nodiscard]] virtual double evaluate(double time, std::size_t component = 0) const = 0;

    /**
     * Evaluate the derivative of the function for the given @p time and @p component.
     */
    [[nodiscard]] virtual double evaluate_derivative(
        double time, std::size_t component = 0) const = 0;
  };

  /**
   * @brief Function based on user-supplied expressions
   *
   * This class supports functions of type \f$ f( t, a_1(t), ..., a_k(t)) \f$, where
   *  \f$ a_1(t), ..., a_k(t) \f$ are time-dependent FunctionVariable objects.
   */
  class SymbolicFunctionOfTime : public FunctionOfTime
  {
   public:
    /**
     * Create a SymbolicFunctionOfTime From a vector of @p expressions and a vector of @p variables.
     * Any time-dependent variables basing on the FunctionVariable must be passed in the @p
     * variables vector.
     */
    SymbolicFunctionOfTime(const std::vector<std::string>& expressions,
        std::vector<Teuchos::RCP<FunctionVariable>> variables);

    [[nodiscard]] double evaluate(double time, std::size_t component = 0) const override;

    [[nodiscard]] double evaluate_derivative(double time, std::size_t component = 0) const override;

   private:
    using ValueType = double;
    using FirstDerivativeType = Sacado::Fad::DFad<double>;

    //! vector of parsed expressions
    std::vector<Teuchos::RCP<Core::UTILS::SymbolicExpression<ValueType>>> expr_;


    //! vector of the function variables and all their definitions
    std::vector<Teuchos::RCP<FunctionVariable>> variables_;
  };

  //! create a vector function of time from multiple expressions
  Teuchos::RCP<FunctionOfTime> try_create_function_of_time(
      const std::vector<Input::LineDefinition>& function_line_defs);

}  // namespace Core::UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
