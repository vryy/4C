/*! \file
\brief Interface for functions of time
\level 0
*/

#ifndef BACI_UTILS_FUNCTION_OF_TIME_H
#define BACI_UTILS_FUNCTION_OF_TIME_H

#include "baci_config.H"

#include "baci_discretization_fem_general_utils_polynomial.H"
#include "baci_io_linedefinition.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_utils_exceptions.H"
#include "baci_utils_function_manager.H"
#include "baci_utils_functionvariables.H"
#include "baci_utils_symbolic_expression.H"

#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

#include <complex>
#include <iostream>
#include <vector>

BACI_NAMESPACE_OPEN

namespace CORE::UTILS
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
    [[nodiscard]] virtual double Evaluate(double time, std::size_t component = 0) const = 0;

    /**
     * Evaluate the derivative of the function for the given @p time and @p component.
     */
    [[nodiscard]] virtual double EvaluateDerivative(
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

    [[nodiscard]] double Evaluate(double time, std::size_t component = 0) const override;

    [[nodiscard]] double EvaluateDerivative(double time, std::size_t component = 0) const override;

   private:
    using ValueType = double;
    using FirstDerivativeType = Sacado::Fad::DFad<double>;

    //! vector of parsed expressions
    std::vector<Teuchos::RCP<CORE::UTILS::SymbolicExpression<ValueType>>> expr_;


    //! vector of the function variables and all their definitions
    std::vector<Teuchos::RCP<FunctionVariable>> variables_;
  };

  //! create a vector function of time from multiple expressions
  Teuchos::RCP<FunctionOfTime> TryCreateFunctionOfTime(
      const std::vector<INPUT::LineDefinition>& function_line_defs);

}  // namespace CORE::UTILS

BACI_NAMESPACE_CLOSE

#endif
