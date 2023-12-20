/*----------------------------------------------------------------------*/
/*! \file

\brief Time dependent variables for function definition

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UTILS_FUNCTIONVARIABLES_H
#define BACI_UTILS_FUNCTIONVARIABLES_H

#include "baci_config.H"

#include "baci_utils_symbolic_expression.H"

#include <Sacado.hpp>

BACI_NAMESPACE_OPEN

// define periodicstruct
struct periodicstruct
{
  bool periodic;
  double t1;
  double t2;
};



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
namespace CORE::UTILS
{
  /// class that defines the main properties of a time dependent variable
  class FunctionVariable
  {
   public:
    /// Base constructor takign the @p name of the FunctionVariable.
    FunctionVariable(std::string name);

    /// Virtual destructor.
    virtual ~FunctionVariable() = default;

    /// evaluate the function at given time
    virtual double Value(const double t) = 0;

    /// evaluate the time derivative at given time
    virtual double TimeDerivativeValue(const double t, const unsigned deg = 1) = 0;

    /// check the inclusion of the considered time
    virtual bool ContainTime(const double t) = 0;

    /// return the name of the variable
    std::string Name() { return name_; }

   private:
    /// name of the variable
    const std::string name_;
  };


  /// variable based on a parsed expression
  class ParsedFunctionVariable : public FunctionVariable
  {
   public:
    ParsedFunctionVariable(std::string name, const std::string& buf);

    /// evaluate the function at given time
    double Value(const double t) override;

    /// evaluate the time derivative at given time
    double TimeDerivativeValue(const double t, const unsigned deg = 1) override;

    /// check the inclusion of the considered time
    bool ContainTime(const double t) override;

   private:
    /// parsed function
    Teuchos::RCP<CORE::UTILS::SymbolicExpression<double>> timefunction_;
  };


  /// variable based on a linear interpolation
  class LinearInterpolationVariable : public FunctionVariable
  {
   public:
    LinearInterpolationVariable(std::string name, std::vector<double> times,
        std::vector<double> values, struct periodicstruct periodicdata);

    /// evaluate the function at given time
    double Value(const double t) override;

    /// templated function to evaluate and to derive from using sacado
    template <typename ScalarT>
    ScalarT Value(const ScalarT& t);

    /// evaluate the time derivative at given time
    double TimeDerivativeValue(const double t, const unsigned deg = 1) override;

    /// check the inclusion of the considered time
    bool ContainTime(const double t) override;

   private:
    /// times for the interpolation
    const std::vector<double> times_;

    /// values for the interpolation
    const std::vector<double> values_;

    /// flag for periodic repetition
    const bool periodic_;

    /// initial time of the periodic repetition
    const double t1_;

    /// final time of the periodic repetition
    const double t2_;
  };


  /// variable based on a set of parsed expressions
  class MultiFunctionVariable : public FunctionVariable
  {
   public:
    MultiFunctionVariable(std::string name, std::vector<double> times,
        std::vector<std::string> description_vec, struct periodicstruct periodicdata);

    /// evaluate the function at given time
    double Value(const double t) override;

    /// evaluate the time derivative at given time
    double TimeDerivativeValue(const double t, const unsigned deg = 1) override;

    /// check the inclusion of the considered time
    bool ContainTime(const double t) override;

   private:
    /// times defining each interval
    const std::vector<double> times_;

    /// vector of parsed functions
    std::vector<Teuchos::RCP<CORE::UTILS::SymbolicExpression<double>>> timefunction_;


    /// flag for periodic repetition
    const bool periodic_;

    /// initial time of the periodic repetition
    const double t1_;

    /// final time of the periodic repetition
    const double t2_;
  };


  /// variable based on a Fourier interpolation
  class FourierInterpolationVariable : public FunctionVariable
  {
   public:
    FourierInterpolationVariable(std::string name, std::vector<double> times,
        std::vector<double> values, struct periodicstruct periodicdata);

    /// evaluate the function at given time
    double Value(const double t) override;

    /// templated function to evaluate and to derive from using sacado
    template <typename ScalarT>
    ScalarT Value(const ScalarT& t);

    /// evaluate the time derivative at given time
    double TimeDerivativeValue(const double t, const unsigned deg = 1) override;

    /// check the inclusion of the considered time
    bool ContainTime(const double t) override;

   private:
    /// times for the interpolation
    const std::vector<double> times_;

    /// values for the interpolation
    const std::vector<double> values_;

    /// flag for periodic repetition
    const bool periodic_;

    /// initial time of the periodic repetition
    const double t1_;

    /// final time of the periodic repetition
    const double t2_;
  };


  /**
   * @brief A FunctionVariable constructed piece-wise from other FunctionVariables.
   *
   * When the function is evaluated with either a call to Value() or TimeDerivativeValue(), the
   * first piece that contains the given time will be used for evaluation, even if multiple pieces
   * would be able to evaluate the time.
   */
  class PiecewiseVariable : public FunctionVariable
  {
   public:
    //! Create a PiecewiseVariable from the given @p pieces.
    PiecewiseVariable(const std::string& name, std::vector<Teuchos::RCP<FunctionVariable>> pieces);

    double Value(double t) override;

    double TimeDerivativeValue(double t, unsigned int deg) override;

    bool ContainTime(double t) override;

   private:
    //! Helper function to access the piece that contains time @p time. Returns the first
    //! FunctionVariable that contains the time @p t and does not check if other FunctionVariables
    //! also contain it.
    FunctionVariable& FindPieceForTime(double t);

    //! Store the pieces that make up the variable in different time intervals.
    std::vector<Teuchos::RCP<FunctionVariable>> pieces_;
  };

}  // namespace CORE::UTILS

BACI_NAMESPACE_CLOSE

#endif
