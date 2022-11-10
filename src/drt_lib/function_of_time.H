/*! \file
\brief Interface for functions of time
\level 0
*/

#ifndef SRC_DRT_LIB_FUNCTION_OF_TIME_H
#define SRC_DRT_LIB_FUNCTION_OF_TIME_H

#include <iostream>
#include <vector>
#include <complex>
#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

#include "linalg_fixedsizematrix.H"
#include "drt_dserror.H"
#include "drt_parser.H"
#include "drt_functionvariables.H"
#include "drt_utils_polynomial.H"
#include "drt_function_manager.H"

namespace DRT::UTILS
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
    std::vector<Teuchos::RCP<DRT::PARSER::Parser<ValueType>>> expr_;

    //! vector of parsed first derivatives
    std::vector<Teuchos::RCP<DRT::PARSER::Parser<FirstDerivativeType>>> dexprdt_;

    //! vector of the function variables and all their definitions
    std::vector<Teuchos::RCP<FunctionVariable>> variables_;
  };

  //! create a vector function of time from multiple expressions
  Teuchos::RCP<FunctionOfTime> TryCreateFunctionOfTime(
      std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> function_lin_defs);

}  // namespace DRT::UTILS

#endif
