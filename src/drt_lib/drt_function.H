/*----------------------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of space- and/or time-dependent functions

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef DRT_FUNCTION_H
#define DRT_FUNCTION_H


#include <iostream>
#include <vector>
#include <complex>
#include <Teuchos_RCP.hpp>

#include "../linalg/linalg_fixedsizematrix.H"
#include "drt_dserror.H"
#include "drt_parser.H"
#include "drt_functionvariables.H"
#include "../drt_fem_general/drt_utils_polynomial.H"

namespace MAT
{
  class Material;
}

namespace Sacado
{
  namespace Fad
  {
    template <typename>
    class DFad;
  }
}  // namespace Sacado

namespace DRT
{
  class Condition;
  class Problem;


  namespace INPUT
  {
    class Lines;
    class DatFileReader;

  }  // namespace INPUT

  namespace UTILS
  {
    /*!

    \brief abstract interface for space- and time-dependent functions.

    Functions can be defined either scalar-valued or vector-valued.

    There exist function definitions using explicitly defined spatial functions
    for test problems (beltrami, kim-moin etc).

    Spatial function based on parsed expressions are provided as well.

    */

    class Function
    {
     public:
      virtual ~Function(){};

      /// Evaluate time and space dependent FUNCT
      virtual double Evaluate(const int index,  ///< component of function
          const double* x,                      ///< coordinates
          double t                              ///< time
          ) = 0;


      /// Evaluate time dependent FUNCT
      double EvaluateTime(double t  ///< time
      )
      {
        double dummy[3] = {0.0};

#ifdef DEBUG
        double dummy2[3] = {123.456};
        if (fabs(Evaluate(0, dummy, t) - Evaluate(0, dummy2, t)) > 1.0e-14)
          dserror("You evaluate a FUNCT solely by time even if it dependents on x,y or z!");
#endif

        return Evaluate(0, dummy, t);
      }


      /// Evaluate first spatial derivatives of time and space dependent FUNCT
      virtual std::vector<double> EvaluateSpatialDerivative(
          const int index,  ///< component of function
          const double* x,  ///< coordinates
          const double t    ///< time
      )
      {
        dserror("The evaluation of the derivative is not implemented for this function");
        std::vector<double> emptyvector;
        return emptyvector;
      };

      /// Evaluate first or second time derivative of time and space dependent FUNCT
      virtual std::vector<double> EvaluateTimeDerivative(
          const int index,    ///< component of function
          const double* x,    ///< coordinates
          const double t,     ///< time
          const unsigned deg  ///< degree of maximal derivative
      )
      {
        dserror("The evaluation of the time derivative is not implemented for this function");
        std::vector<double> emptyvector;
        return emptyvector;
      };


      /// Evaluate first or second time derivative of time dependent FUNCT
      std::vector<double> EvaluateTimeDerivative(const double t,  ///< time
          const unsigned deg                                      ///< degree of maximal derivative
      )
      {
        double dummy[3] = {0.0};

#ifdef DEBUG
        std::vector<double> res1 = EvaluateTimeDerivative(0, dummy, t, deg);
        double dummy2[3] = {123.456};
        std::vector<double> res2 = EvaluateTimeDerivative(0, dummy2, t, deg);
        for (unsigned int i = 0; i < deg; i++)
        {
          if (fabs(res1[deg] - res2[deg]) > 1.0e-14)
            dserror(
                "You evaluate the time derivative of a FUNCT solely by time even if it dependents "
                "on x,y or z!");
        }
#endif

        return EvaluateTimeDerivative(0, dummy, t, deg);
      };


      /// Return number of components of FUNCT
      virtual int NumberComponents() { return 0; }
    };



    class ExprFunction : public Function
    {
     public:
      ExprFunction();

      virtual ~ExprFunction();

      /// add an expression to the function
      virtual void AddExpr(
          std::string buf, std::vector<std::vector<Teuchos::RCP<FunctionVariable>>> variables);

      /// evaluate the function at given space and time
      virtual double Evaluate(const int index, const double* x, double t);

      /// evaluate the first derivative of the function at given space and time
      virtual std::vector<double> EvaluateSpatialDerivative(
          const int index, const double* x, const double t);

      /// evaluate the 1st or 2nd time derivative of the function at given space and time
      virtual std::vector<double> EvaluateTimeDerivative(
          const int index, const double* x, const double t, const unsigned deg);

      /// return the number of components
      virtual int NumberComponents() { return (expr_.size()); }

      /// check if it is a variable
      bool IsVariable(const int index, const std::string& varname) const;

      /// add a variable to the function
      void AddVariable(const int index, const std::string& varname, double varvalue);

      /// pars the expressions
      void ParseExpressions();

     protected:
      /// dimension of the problem
      int dim_;

      /// vector of parsed expressions
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>> expr_;

      /// vector of parsed 1st derivatives
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<double>>>> exprd_;

      /// vector of parsed 2nd derivatives
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>>>
          exprdd_;

      /// vector of the function variables and all their definitions
      std::vector<std::vector<Teuchos::RCP<FunctionVariable>>> variables_;

     private:
      bool isparsed_;
    };



    /// class that manages all defined spatial functions
    class FunctionManager
    {
     public:
      FunctionManager() {}

      Teuchos::RCP<DRT::INPUT::Lines> ValidFunctionLines();

      /// read the baci input file and set up all functions
      void ReadInput(DRT::INPUT::DatFileReader& reader);

      /// return a specified function
      Function& Funct(int num);

     private:
      /// no copying
      FunctionManager(const FunctionManager&);

      /// no copying
      FunctionManager& operator=(const FunctionManager&);

      /// all functions defined. The functions may be vector valued!
      std::vector<Teuchos::RCP<Function>> functions_;
    };


    /// spatial function based on parsed expression
    class VariableExprFunction : public ExprFunction
    {
     public:
      /*!

      \brief Default constructor creating empty object. Expressions are
             added with add function

      */
      VariableExprFunction();

      /*!

      \brief clean up parse tree

      */
      virtual ~VariableExprFunction();

      /*!

      \brief add expression to an existing ExprFunction in order to extend
             it to a vector-valued spatial function.

             Every call to AddExpr adds one more component to the
             vector-valued function.

      \param buf (i) (c-string) expression to be parsed during evaluation of this component
      \param constants   (i) vector containing constants defined by a pair (name and value)

      */
      virtual void AddExpr(std::string buf, std::vector<std::pair<std::string, double>> constants);

      /*!

      \brief evaluate function at given time and position

      \param index (i) For vector-valued functions, index defines the
                       function-component which should be evaluated
                       For scalar functions, index is always set to 0
      \param t     (i) The point in time in which the
                       function will be evaluated
      \param x     (i) The point in 3-dimensional space in which the
                       function will be evaluated

      */
      virtual double Evaluate(const int index, const double* x, const double t);

      /*!

      \brief evaluate function for a given set of variables

      \param index    (i)  For vector-valued functions, index defines the
                           function-component which should be evaluated
      \param variables (i) A vector containing a pair (variablename, value) for each variable

      */
      virtual double Evaluate(
          const int index, const std::vector<std::pair<std::string, double>>& variables);

      /*!

      \brief evaluate function for a given set of variables

      \param index    (i)  For vector-valued functions, index defines the
                           function-component which should be evaluated
      \param variables (i) A vector containing a pair (variablename, value) for each variable
      \param constants (i) A vector containing a pair (variablename, value) for each constant

      */
      virtual double Evaluate(const int index,
          const std::vector<std::pair<std::string, double>>& variables,
          const std::vector<std::pair<std::string, double>>& constants);

      /*!

      \brief evaluate derivatives function at given position in space

      \param index (i) For vector-valued functions, index defines the
                       function-component which should be evaluated
                       For scalar functions, index is always set to 0
      \param x     (i) The point in 3-dimensional space in which the
                       function will be evaluated
      \param t     (i) Absolute time in which the
                       function will be evaluated

      */
      virtual std::vector<double> EvaluateSpatialDerivative(
          const int index, const double* x, const double t);

      /*!

      \brief  at given position in space for a given set of variables

      \param index    (i)  For vector-valued functions, index defines the
                           function-component which should be evaluated
      \param variables (i) A vector containing a pair (variablename, value) for each variable

      */
      virtual std::vector<double> EvaluateDerivative(
          const int index, const std::vector<std::pair<std::string, double>>& variables);

      /*!

      \brief  at given position in space for a given set of variables and constants

      \param index    (i)  For vector-valued functions, index defines the
                           function-component which should be evaluated
      \param variables (i) A vector containing a pair (variablename, value) for each variable
      \param constants (i) A vector containing a pair (variablename, value) for each constant

      */
      virtual std::vector<double> EvaluateDerivative(int index,
          const std::vector<std::pair<std::string, double>>& variables,
          const std::vector<std::pair<std::string, double>>& constants);

      /*!
      \brief  check if a given name is a variable

      \param index    (i)  For vector-valued functions, index defines the
                           function-component which should be evaluated
      \param varname  (i)  Variable name to be checked

      */
      virtual bool IsVariable(const int index, const std::string& varname) const;

      /*!
      \brief  add a variable to the parser object

      \param index    (i)  For vector-valued functions, index defines the
                           function-component which should be evaluated
      \param varname  (i)  Variable name

      \param varvalue (i)  Value of the variable

      */
      virtual void AddVariable(const int index, const std::string& varname, double varvalue);

      /*!
      \brief  parse the underlying expression
      */
      virtual void ParseExpressions();

     private:
      //! flag indicating whether expression was parsed
      bool isparsed_;
    };

  }  // namespace UTILS

}  // namespace DRT


void PrintFunctionDatHeader();


#endif
