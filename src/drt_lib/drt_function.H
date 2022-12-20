/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluating of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef DRT_FUNCTION_H
#define DRT_FUNCTION_H


#include <iostream>
#include <vector>
#include <complex>
#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

#include "linalg_fixedsizematrix.H"
#include "drt_dserror.H"
#include "drt_functionvariables.H"
#include "drt_utils_polynomial.H"
#include "drt_function_manager.H"


namespace DRT
{

  namespace UTILS
  {

    template <class T>
    class SymbolicExpression;

    class FunctionManager;

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
      virtual double Evaluate(const double* x, double t, std::size_t component) = 0;

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
          const double* x, double t, std::size_t component)
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
          const double* x, double t, unsigned deg, std::size_t component)
      {
        dserror("The evaluation of the time derivative is not implemented for this function");
        std::vector<double> emptyvector;
        return emptyvector;
      };

      /// Return number of components of function
      virtual std::size_t NumberComponents() { return 0; }
    };


    /**
     * @brief Function based on user-supplied expressions
     *
     * This class supports functions of type \f$ f(\mathbf{x}, t, a_1(t), ..., a_k(t)) \f$, where
     * \f$ \mathbf{x} \f$ is the spatial coordinate with `dim` components and where \f$ a_1(t), ...,
     * a_k(t) \f$ are time-dependent FunctionVariable objects.
     */
    template <int dim>
    class ExprFunction : public FunctionOfSpaceTime
    {
     public:
      /**
       * Create an ExprFunction. Each entry in @p expressions corresponds to one component of the
       * function, thus the resulting function will have `expressions.size()` components. Any
       * time-dependent variables that appear in the expressions must be passed in the @p variables
       * vector.
       */
      ExprFunction(const std::vector<std::string>& expressions,
          std::vector<Teuchos::RCP<FunctionVariable>> variables);

      double Evaluate(const double* x, double t, std::size_t component) override;

      std::vector<double> EvaluateSpatialDerivative(
          const double* x, double t, std::size_t component) override;

      std::vector<double> EvaluateTimeDerivative(
          const double* x, double t, unsigned deg, std::size_t component) override;

      std::size_t NumberComponents() override { return (expr_.size()); }

     private:
      using ValueType = double;
      using SecondDerivativeType = Sacado::Fad::DFad<Sacado::Fad::DFad<ValueType>>;

      /// vector of parsed expressions
      std::vector<Teuchos::RCP<DRT::UTILS::SymbolicExpression<ValueType>>> expr_;

      /// vector of the function variables and all their definitions
      std::vector<Teuchos::RCP<FunctionVariable>> variables_;
    };


    /// vector-valued spatial function based on parsed expression
    template <int dim>
    class VariableExprFunction : public FunctionOfSpaceTime
    {
     public:
      VariableExprFunction(
          const std::string& component, std::vector<std::pair<std::string, double>> constants);

      VariableExprFunction() = default;

      double Evaluate(const double* x, double t, std::size_t component) override;



      /*!
       * \brief evaluate function for a given set of variables
       *
       * \param index (i) For vector-valued functions, index defines the function-component which
       * should be evaluated
       * \param variables (i) A vector containing a pair (variablename, value)
       * for each variable
       */
      virtual double Evaluate(
          const int index, const std::vector<std::pair<std::string, double>>& variables);

      /*!
       * \brief evaluate function for a given set of variables
       *
       * \param index (i) For vector-valued functions, index defines the function-component which
       * should be evaluated
       * \param variables (i) A vector containing a pair (variablename, value)
       * for each variable
       * \param constants (i) A vector containing a pair (variablename, value) for
       * each constant
       */
      virtual double Evaluate(const int index,
          const std::vector<std::pair<std::string, double>>& variables,
          const std::vector<std::pair<std::string, double>>& constants);

      std::vector<double> EvaluateSpatialDerivative(
          const double* x, double t, std::size_t component) override;

      /*!
       *\brief evaluate the derivative at given position in space for a given set of variables
       *
       * \param index (i) For vector-valued functions, index defines the function-component which
       * should be evaluated
       * \param variables (i) A vector containing a pair (variablename, value)
       * for each variable
       */
      virtual std::vector<double> EvaluateDerivative(
          const int index, const std::vector<std::pair<std::string, double>>& variables);

      /*!
       * \brief evaluate the derivative at given position in space for a given set of variables and
       * constants
       *
       * \param index (i) For vector-valued functions, index defines the function-component which
       * should be evaluated
       * \param variables (i) A vector containing a pair (variablename, value)
       * for each variable
       * \param constants (i) A vector containing a pair (variablename, value) for
       * each constant
       */
      virtual std::vector<double> EvaluateDerivative(int index,
          const std::vector<std::pair<std::string, double>>& variables,
          const std::vector<std::pair<std::string, double>>& constants);


      /// return the number of components
      std::size_t NumberComponents() override { return (expr_.size()); }

     private:
      using ValueType = double;
      using FirstDerivativeType = Sacado::Fad::DFad<ValueType>;

      /// dimension of the problem
      int dim_;

      /// vector of parsed expressions
      std::vector<Teuchos::RCP<DRT::UTILS::SymbolicExpression<ValueType>>> expr_;

      /// vector of the function variables and all their definitions
      std::vector<std::vector<Teuchos::RCP<FunctionVariable>>> variables_;

     private:
      // constants from input
      std::vector<std::pair<std::string, ValueType>> constants_from_input_;

      //! throw an error if a constant given in the input file is a primary variables
      template <typename T>
      void CheckValidInput(const std::map<std::string, T>& variable_values) const;
    };

    /// try to create VariableExprFunction from a given line definition
    template <int dim>
    Teuchos::RCP<FunctionOfSpaceTime> TryCreateVariableExprFunction(
        Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def,
        DRT::UTILS::FunctionManager& manager, const int index_current_funct_in_manager);

    /// create a vector function from multiple expressions
    template <int dim>
    Teuchos::RCP<FunctionOfSpaceTime> TryCreateExprFunction(
        std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> functions_lin_defs);
  }  // namespace UTILS

}  // namespace DRT

#endif
