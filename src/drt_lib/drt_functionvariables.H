/*----------------------------------------------------------------------*/
/*!
\file drt_functionvariables.H

\brief Time dependent variables for function definition

<pre>
\level 1

\maintainer Andrea La Spina
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef DRT_FUNCTIONVARIABLES_H
#define DRT_FUNCTIONVARIABLES_H

// define periodicstruct
struct periodicstruct
{
  bool periodic;
  double t1;
  double t2;
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
namespace Sacado
{
namespace Fad
{
template <typename> class DFad;
}
}

namespace DRT {

namespace PARSER
{
template <typename> class Parser;
}

namespace UTILS {

/// class that defines the main properties of a time dependent variable
class FunctionVariable
{
public:
  /// evaluate the function at given time
  virtual double Value(const double t=0) = 0;

  /// evaluate the time derivative at given time
  virtual double TimeDerivativeValue(const double t=0, const unsigned deg=1) = 0;

  /// check the inclusion of the considered time
  virtual bool ContainTime (const double t=0) { return false; };

  /// return the name of the variable
  std::string Name() { return name_; }

protected:
  /// name of the variable
  std::string name_;

private:
};


/// variable based on a parsed expression
class ParsedFunctionVariable : public FunctionVariable
{
public:
  ParsedFunctionVariable(std::string name, std::string buf);

  ~ParsedFunctionVariable();

  /// evaluate the function at given time
  virtual double Value(const double t);

  /// evaluate the time derivative at given time
  virtual double TimeDerivativeValue(const double t, const unsigned deg);

  /// check the inclusion of the considered time
  virtual bool ContainTime (const double t);

private:
  /// parsed function
  Teuchos::RCP<DRT::PARSER::Parser<double> > timefunction_;

  // parsed derivative
  Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > > > timederivative_;
};


/// variable based on a linear interpolation
class LinearInterpolationVariable : public FunctionVariable
{
public:
  LinearInterpolationVariable(std::string name, std::vector<double> times, std::vector<double> values, struct periodicstruct periodicdata);

  ~LinearInterpolationVariable();

  /// evaluate an equivalent time for a periodic repetition
  double EquivalentTime(const double t);

  /// evaluate the function at given time
  virtual double Value(const double t);

  /// evaluate the time derivative at given time
  virtual double TimeDerivativeValue(const double t, const unsigned deg);

  /// check the inclusion of the considered time
  virtual bool ContainTime (const double t);

private:
  /// times for the interpolation
  std::vector<double> times_;

  /// values for the interpolation
  std::vector<double> values_;

  /// flag for periodic repetition
  bool periodic_;

  /// initial time of the periodic repetition
  double t1_;

  /// final time of the periodic repetition
  double t2_;
};


/// variable based on a set of parsed expressions
class MultiFunctionVariable : public FunctionVariable
{
public:
  MultiFunctionVariable(std::string name, std::vector<double> times, std::vector<std::string> description_vec, struct periodicstruct periodicdata);

  ~MultiFunctionVariable();

  /// evaluate an equivalent time for a periodic repetition
  double EquivalentTime(const double t);

  /// evaluate the function at given time
  virtual double Value(const double t);

  /// evaluate the time derivative at given time
  virtual double TimeDerivativeValue(const double t, const unsigned deg);

  /// check the inclusion of the considered time
  virtual bool ContainTime (const double t);

private:
  /// times defining each interval
  std::vector<double> times_;

  ///vector of parsed functions
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<double> > > timefunction_;

  /// vector of parsed derivatives
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > > > > timederivative_;

  /// flag for periodic repetition
  bool periodic_;

  /// initial time of the periodic repetition
  double t1_;

  /// final time of the periodic repetition
  double t2_;
};


/// variable based on a Fourier interpolation
class FourierInterpolationVariable : public FunctionVariable
{
public:
  FourierInterpolationVariable(std::string name, std::vector<double> times, std::vector<double> values, struct periodicstruct periodicdata);

  ~FourierInterpolationVariable();

  /// evaluate an equivalent time for a periodic repetition
  double EquivalentTime(const double t);

  /// evaluate the function at given time
  virtual double Value(const double t);
  template <typename ScalarT>
  ScalarT Value(const ScalarT& t);

  /// evaluate the time derivative at given time
  virtual double TimeDerivativeValue(const double t, const unsigned deg);

  /// check the inclusion of the considered time
  virtual bool ContainTime (const double t);

private:
  /// times for the interpolation
  std::vector<double> times_;

  /// values for the interpolation
  std::vector<double> values_;

  /// flag for periodic repetition
  bool periodic_;

  /// initial time of the periodic repetition
  double t1_;

  /// final time of the periodic repetition
  double t2_;
};


} // namespace DRT
} // namespace UTILS

#endif
