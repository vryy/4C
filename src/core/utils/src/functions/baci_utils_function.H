/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluating of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UTILS_FUNCTION_H
#define BACI_UTILS_FUNCTION_H


#include "baci_discretization_fem_general_utils_polynomial.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_utils_exceptions.H"
#include "baci_utils_functionvariables.H"

#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

#include <vector>

BACI_NAMESPACE_OPEN

namespace DRT::INPUT
{
  class LineDefinition;
}

namespace CORE::UTILS
{

  template <class T>
  class SymbolicExpression;

  /*!
   * \brief interface for space- and time-dependent functions.
   *
   * Functions can be defined either scalar-valued or vector-valued.
   */
  class FunctionOfSpaceTime
  {
   public:
    //! Virtual destructor.
    virtual ~FunctionOfSpaceTime() = default;

    /*!
     * @brief Evaluation of time and space dependent function
     *
     * Evaluate the specified component of the function at the specified position and point in
     * time.
     *
     * @param x  (i) The point in 3-dimensional space in which the function will be evaluated
     * @param t  (i) The point in time in which the function will be evaluated
     * @param component (i) For vector-valued functions, index defines the function-component
     *                      which should be evaluated
     * @return function value
     */
    virtual double Evaluate(const double* x, double t, std::size_t component) const = 0;

    /*!
     * \brief Evaluation of first spatial derivative of time and space dependent function
     *
     * \param x  (i) The point in 3-dimensional space in which the function will be evaluated
     * \param t  (i) The point in time in which the function will be evaluated
     * @param component (i) For vector-valued functions, index defines the function-component
     *                      which should be evaluated
     * \return first spatial derivative of function
     */
    virtual std::vector<double> EvaluateSpatialDerivative(
        const double* x, double t, std::size_t component) const
    {
      dserror("The evaluation of the derivative is not implemented for this function");
      std::vector<double> emptyvector;
      return emptyvector;
    };

    /*!
     * \brief Evaluation of time derivatives and value of the time and space dependent function
     *
     * Evaluate the specified component of the function at the specified position and point in
     * time and calculate the time derivative(s) up to degree @p deg.
     *
     * @param x  (i) The point in 3-dimensional space in which the function will be evaluated
     * @param t  (i) The point in time in which the function will be evaluated
     * @param deg   (i) maximum time derivative degree
     * @param component (i) For vector-valued functions, index defines the function-component
     *                      which should be evaluated
     * @return vector containing value and time derivative(s)
     */
    virtual std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const
    {
      dserror("The evaluation of the time derivative is not implemented for this function");
      std::vector<double> emptyvector;
      return emptyvector;
    };

    /// Return number of components of function
    virtual std::size_t NumberComponents() const { return 0; }
  };


  /**
   * @brief Function based on user-supplied expressions
   *
   * This class supports functions of type \f$ f(\mathbf{x}, t, a_1(t), ..., a_k(t)) \f$, where
   * \f$ \mathbf{x} \f$ is the spatial coordinate with `dim` components and where \f$ a_1(t), ...,
   * a_k(t) \f$ are time-dependent FunctionVariable objects.
   */
  template <int dim>
  class SymbolicFunctionOfSpaceTime : public FunctionOfSpaceTime
  {
   public:
    /**
     * Create an SymbolicFunctionOfSpaceTime. Each entry in @p expressions corresponds to one
     * component of the function, thus the resulting function will have `expressions.size()`
     * components. Any time-dependent variables that appear in the expressions must be passed in
     * the @p variables vector.
     */
    SymbolicFunctionOfSpaceTime(const std::vector<std::string>& expressions,
        std::vector<Teuchos::RCP<FunctionVariable>> variables);

    double Evaluate(const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateSpatialDerivative(
        const double* x, double t, std::size_t component) const override;

    std::vector<double> EvaluateTimeDerivative(
        const double* x, double t, unsigned deg, std::size_t component) const override;

    std::size_t NumberComponents() const override { return (expr_.size()); }

   private:
    using ValueType = double;
    using SecondDerivativeType = Sacado::Fad::DFad<Sacado::Fad::DFad<ValueType>>;

    /// vector of parsed expressions
    std::vector<Teuchos::RCP<CORE::UTILS::SymbolicExpression<ValueType>>> expr_;

    /// vector of the function variables and all their definitions
    std::vector<Teuchos::RCP<FunctionVariable>> variables_;
  };


  /*!
   * \brief Interface for mathematical functions and with arbitrary arguments
   *
   * Functions that derive from this interface are free to take arbitrary arguments for
   * evaluation. It is rather obvious that such an interface would encompass all specialized
   * interfaces for functions. Indeed, any function could be implemented under this interface.
   * However, a lot of our functions have a clear interface, e.g. a function depending only on
   * space and time. There are specialized interfaces for those functions and the present
   * interface should not be misused for these cases. This interface is best used for convenience
   * when prototyping new equations. If a new function type emerges during this process it is
   * likely a good idea to create a new interface for it. Since evaluation of a FunctionOfAnything
   * must always go through another level of indirection due to the data structure for the passed
   * variables, there are two pitfalls: one can either pass variables (from the outside) that are
   * not needed or try to access variables that are not passed (inside the function). Derived
   * classes need to deal with this problem on their own. An important derived class of this
   * interface is the SymbolicFunctionOfAnything.
   *
   * This class evaluates and forms the derivative of arbitrary functions for a given set of
   * variables and constants.
   *
   * The functions can be defined either scalar-valued or vector-valued.
   *
   */
  class FunctionOfAnything
  {
   public:
    //! Virtual destructor.
    virtual ~FunctionOfAnything() = default;

    /*!
     * \brief evaluate function for a given set of variables and constants
     *
     * \note there is no distinction between the input arguments variables and constants. So for
     * the function evaluation it makes no difference whether all necessary variables and
     * constants are passed in a single vector together with an empty vector or the variables and
     * constants are passed separately in individual vectors.
     *
     * \param variables (i) A vector containing a pair (variablename, value) for
     * each variable
     * \param constants (i) A vector containing a pair (variablename, value) for each
     * constant
     * \param component (i) For vector-valued functions, component defines the function-component
     * which should be evaluated
     */
    virtual double Evaluate(const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants,
        const std::size_t component) const = 0;

    /*!
     * \brief evaluates the derivative of a function with respect to the given variables
     *
     * \param variables (i) A vector containing a pair (variablename, value) for
     * each variable. The derivative of the function is build with respect to all these variables.
     * \param constants (i) A vector containing a pair (variablename, value) for each
     * constant.
     * \param component (i) For vector-valued functions, component defines the function-component
     * which should be evaluated
     */
    virtual std::vector<double> EvaluateDerivative(
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants,
        const std::size_t component) const = 0;

    //! Return number of components of function
    virtual std::size_t NumberComponents() const = 0;
  };


  /**
   * @brief Function to evaluate and form the derivative of user defined
   * symbolic expressions.
   *
   * The expression must only contain supported functions, literals and operators, as well as
   * arbitrary number of variables and constants. See the class documentation of
   * SymbolicExpression for more details.
   *
   * It is possible to predefine values of constants in the input file.
   */
  template <int dim>
  class SymbolicFunctionOfAnything : public FunctionOfAnything
  {
   public:
    SymbolicFunctionOfAnything(
        const std::string& component, std::vector<std::pair<std::string, double>> constants);


    double Evaluate(const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants,
        const std::size_t component) const override;


    std::vector<double> EvaluateDerivative(
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants,
        const std::size_t component) const override;

    /// return the number of components
    std::size_t NumberComponents() const override { return (expr_.size()); }

   private:
    using ValueType = double;

    //! vector of parsed expressions
    std::vector<Teuchos::RCP<CORE::UTILS::SymbolicExpression<ValueType>>> expr_;

    //! vector of the function variables and all their definitions
    std::vector<std::vector<Teuchos::RCP<FunctionVariable>>> variables_;

   private:
    //! constants from input
    std::vector<std::pair<std::string, ValueType>> constants_from_input_;
  };

  /// try to create SymbolicFunctionOfAnything from a given line definition
  template <int dim>
  Teuchos::RCP<FunctionOfAnything> TryCreateSymbolicFunctionOfAnything(
      const std::vector<BACI::DRT::INPUT::LineDefinition>& function_line_defs);

  /// create a vector function from multiple expressions
  template <int dim>
  Teuchos::RCP<FunctionOfSpaceTime> TryCreateSymbolicFunctionOfSpaceTime(
      const std::vector<BACI::DRT::INPUT::LineDefinition>& function_line_defs);
}  // namespace CORE::UTILS



BACI_NAMESPACE_CLOSE

#endif
