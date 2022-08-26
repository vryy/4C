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

#include "../linalg/linalg_fixedsizematrix.H"
#include "drt_dserror.H"
#include "drt_parser.H"
#include "drt_functionvariables.H"
#include "../drt_fem_general/drt_utils_polynomial.H"
#include "drt_function_manager.H"
#include "function_interface_temporary.h"


namespace DRT
{
  namespace UTILS
  {
    // forward declaration
    class FunctionManager;

    /*!
     * \brief abstract interface for space- and time-dependent functions.
     *
     * Functions can be defined either scalar-valued or vector-valued.
     *
     * There exist function definitions using explicitly defined spatial functions
     * for test problems (beltrami, kim-moin etc).
     *
     * Spatial function based on parsed expressions are provided as well.
     */
    class FunctionOfSpaceTime : public TemporaryFunctionInterface
    {
     public:
      //! Virtual destructor.
      virtual ~FunctionOfSpaceTime() = default;

      /*!
       * @brief Evaluation of time and space dependent FUNCT
       *
       * Evaluate the specified component of the function at the specified position and point in
       * time.
       *
       * @param index (i) For vector-valued functions, index defines the function-component which
       *                  should be evaluated
       * @param x  (i) The point in 3-dimensional space in which the function will be evaluated
       * @param t  (i) The point in time in which the function will be evaluated
       * @return function value
       */
      virtual double Evaluate(const int index, const double* x, double t) = 0;

      /*!
       * @brief Evaluation of time dependent FUNCT
       *
       * @param t  (i) The point in time in which the function will be evaluated
       */
      double EvaluateTime(double t)
      {
        std::array<double, 3> dummy = {0.0};

#ifdef DEBUG
        std::array<double, 3> dummy2 = {123.456};
        if (std::abs(Evaluate(0, dummy.data(), t) - Evaluate(0, dummy2.data(), t)) > 1.0e-14)
          dserror("You evaluate a FUNCT solely by time even if it dependents on x,y or z!");
#endif

        return Evaluate(0, dummy.data(), t);
      }

      /*!
       * \brief Evaluation of first spatial derivatives of time and space dependent FUNCT
       *
       * \param index (i) For vector-valued functions, index defines the function-component which
       *                  should be evaluated. For scalar functions, index is always set to 0
       * \param x  (i) The point in 3-dimensional space in which the function will be evaluated
       * \param t  (i) Absolute time in which the function will be evaluated
       * \return first spatial derivative of FUNCT
       */
      virtual std::vector<double> EvaluateSpatialDerivative(
          const int index, const double* x, const double t)
      {
        dserror("The evaluation of the derivative is not implemented for this function");
        std::vector<double> emptyvector;
        return emptyvector;
      };

      /*!
       * \brief Evaluation of time derivatives and value of the time and space dependent FUNCT
       *
       * Evaluate the specified component of the function at the specified position and point in
       * time and calculate the time derivative(s).
       *
       * @param index (i) For vector-valued functions, index defines the function-component which
       *                  should be evaluated
       * @param x  (i) The point in 3-dimensional space in which the function will be evaluated
       * @param t  (i) The point in time in which the function will be evaluated
       * @param deg   (i) maximum time derivative degree
       * @return vector containing value and time derivative(s)
       */
      virtual std::vector<double> EvaluateTimeDerivative(
          const int index, const double* x, const double t, const unsigned deg)
      {
        dserror("The evaluation of the time derivative is not implemented for this function");
        std::vector<double> emptyvector;
        return emptyvector;
      };

      /*!
       * \brief Evaluation of time derivatives of time dependent FUNCT
       *
       * @param t  (i) The point in time in which the function will be evaluated
       * @param deg   (i) maximum time derivative degree
       * @return vector containing value and time derivative(s)
       */
      std::vector<double> EvaluateTimeDerivative(const double t, const unsigned deg)
      {
        std::array<double, 3> dummy = {0.0};

#ifdef DEBUG
        std::vector<double> res1 = EvaluateTimeDerivative(0, dummy.data(), t, deg);
        std::array<double, 3> dummy2 = {123.456};
        std::vector<double> res2 = EvaluateTimeDerivative(0, dummy2.data(), t, deg);
        for (unsigned int i = 0; i < deg; i++)
        {
          if (std::abs(res1[deg] - res2[deg]) > 1.0e-14)
          {
            dserror(
                "You evaluate the time derivative of a FUNCT solely by time even if it dependents "
                "on x,y or z!");
          }
        }
#endif

        return EvaluateTimeDerivative(0, dummy.data(), t, deg);
      };

      /// Return number of components of FUNCT
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

      double Evaluate(const int index, const double* x, double t) override;

      std::vector<double> EvaluateSpatialDerivative(
          const int index, const double* x, const double t) override;

      std::vector<double> EvaluateTimeDerivative(
          const int index, const double* x, const double t, const unsigned deg) override;

      /// return the number of components
      std::size_t NumberComponents() override { return (expr_.size()); }

     private:
      using ValueType = double;
      using SecondDerivativeType = Sacado::Fad::DFad<Sacado::Fad::DFad<double>>;

      /// vector of parsed expressions
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<ValueType>>> expr_;

      /// vector of parsed 2nd derivatives
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<SecondDerivativeType>>> exprdd_;

      /// vector of the function variables and all their definitions
      std::vector<Teuchos::RCP<FunctionVariable>> variables_;
    };


    /// vector-valued spatial function based on parsed expression
    template <int dim>
    class VariableExprFunction : public FunctionOfSpaceTime
    {
     public:
      /*!
       * \brief add expression to an existing ExprFunction in order to extend it to a vector-valued
       * spatial function.
       *
       * Every call to AddExpr adds one more component to the vector-valued function.
       * \param buf (i) (c-string) expression to be parsed during evaluation of this component
       * \param constants   (i) vector containing constants defined by a pair (name and value)
       */
      virtual void AddExpr(
          const std::string& buf, const std::vector<std::pair<std::string, double>>& constants);

      double Evaluate(const int index, const double* x, const double t) override;

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
          const int index, const double* x, const double t) override;

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

      /*!
       * \brief check if a given name is a variable
       *
       * \param index (i) For vector-valued functions, index defines the function-component
       *                      which should be evaluated
       * \param varname  (i) Variable name to be checked
       */
      virtual bool IsVariable(const int index, const std::string& varname) const;

      /*!
       * \brief add a variable to the function
       *
       * \param index (i) For vector-valued functions, index defines the function-component
       *                      which should be evaluated
       * \param varname  (i) Variable name
       * \param varvalue (i) Value of the variable
       */
      virtual void AddVariable(const int index, const std::string& varname, double varvalue);

      /*!
       * \brief parse the underlying expression
       */
      virtual void ParseExpressions();

      /// return the number of components
      std::size_t NumberComponents() override { return (expr_.size()); }

     private:
      /// flag indicating whether expression was parsed
      bool isparsed_{false};

      /// dimension of the problem
      int dim_;

      /// vector of parsed expressions
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>> expr_;

      /// vector of parsed 1st derivatives
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<double>>>> exprd_;

      /// vector of the function variables and all their definitions
      std::vector<std::vector<Teuchos::RCP<FunctionVariable>>> variables_;
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
