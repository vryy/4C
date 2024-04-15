/*----------------------------------------------------------------------*/
/*! \file

\brief Time dependent variables for function definition

\level 0


*/
/*----------------------------------------------------------------------*/

#include "baci_utils_functionvariables.hpp"

#include "baci_utils_symbolic_expression.hpp"

#include <Sacado.hpp>

#include <utility>

FOUR_C_NAMESPACE_OPEN


CORE::UTILS::FunctionVariable::FunctionVariable(std::string name) : name_(std::move(name)) {}


CORE::UTILS::ParsedFunctionVariable::ParsedFunctionVariable(
    std::string name, const std::string& buf)
    : FunctionVariable(std::move(name)),
      timefunction_(Teuchos::rcp(new CORE::UTILS::SymbolicExpression<double>(buf)))

{
}


double CORE::UTILS::ParsedFunctionVariable::Value(const double t)
{
  // evaluate the value of the function
  double value = timefunction_->Value({{"t", t}});

  return value;
}


double CORE::UTILS::ParsedFunctionVariable::TimeDerivativeValue(const double t, const unsigned deg)
{
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> tfad(1, 0, t);
  tfad.val() = Sacado::Fad::DFad<double>(1, 0, t);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> vfad;

  vfad = timefunction_->SecondDerivative({{"t", tfad}}, {});

  if (deg == 0)
  {
    const double v = Value(t);
    return v;
  }
  if (deg == 1)
  {
    const double dv_dt = vfad.dx(0).val();
    return dv_dt;
  }
  else if (deg == 2)
  {
    const double d2v_dt2 = vfad.dx(0).dx(0);
    return d2v_dt2;
  }
  else
  {
    dserror("Higher than second derivative is not implemented!");
    return 0.0;
  }
}


bool CORE::UTILS::ParsedFunctionVariable::ContainTime(const double t) { return true; }


CORE::UTILS::LinearInterpolationVariable::LinearInterpolationVariable(std::string name,
    std::vector<double> times, std::vector<double> values, struct Periodicstruct periodicdata)
    : FunctionVariable(std::move(name)),
      times_(std::move(times)),
      values_(std::move(values)),
      periodic_(periodicdata.periodic),
      t1_(periodicdata.t1),
      t2_(periodicdata.t2)
{
}


double CORE::UTILS::LinearInterpolationVariable::Value(const double t)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t >= t1_ - 1.0e-14 and t <= t2_ + 1.0e-14)
  {
    t_equivalent = fmod(t + 1.0e-14, times_[times_.size() - 1] - times_[0]) - 1.0e-14;
  }

  const auto value = Value<double>(t_equivalent);

  return value;
}

template <typename ScalarT>
ScalarT CORE::UTILS::LinearInterpolationVariable::Value(const ScalarT& t)
{
  ScalarT value = 0.0;

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t < times_[0] - 1.0e-14)
  {
    dserror("t_equivalent is smaller than the starting time");
    return 0.0;
  }
  else if (t <= times_[0] + 1.0e-14)
  {
    index = 1;
  }
  else
  {
    while (t_temp <= t - 1.0e-14)
    {
      index++;
      t_temp = times_[index];

      if (index == times_.size())
      {
        dserror("t_equivalent is bigger than the ending time");
        return 0.0;
      }
    }
  }

  // evaluate the value of the function
  value = values_[index - 1] + (values_[index] - values_[index - 1]) /
                                   (times_[index] - times_[index - 1]) * (t - times_[index - 1]);

  return value;
}


double CORE::UTILS::LinearInterpolationVariable::TimeDerivativeValue(
    const double t, const unsigned deg)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t >= t1_ - 1.0e-14 and t <= t2_ + 1.0e-14)
  {
    t_equivalent = fmod(t + 1.0e-14, times_[times_.size() - 1] - times_[0]) - 1.0e-14;
  }

  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> tfad(1, 0, t_equivalent);
  tfad.val() = Sacado::Fad::DFad<double>(1, 0, t_equivalent);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> vfad =
      Value<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>(tfad);

  if (deg == 0)
  {
    const double v = Value(t);
    return v;
  }
  if (deg == 1)
  {
    const double dv_dt = vfad.dx(0).val();
    return dv_dt;
  }
  else if (deg == 2)
  {
    //    const double d2v_dt2 = vfad.dx(0).dx(0);
    return 0.0;
  }
  else
  {
    dserror("Higher than second derivative is not implemented!");
    return 0.0;
  }
}


bool CORE::UTILS::LinearInterpolationVariable::ContainTime(const double t)
{
  /// check the inclusion of the considered time
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t >= t1_ - 1.0e-14 and t <= t2_ + 1.0e-14)
  {
    t_equivalent = fmod(t + 1.0e-14, times_[times_.size() - 1] - times_[0]) - 1.0e-14;
  }

  if (t_equivalent >= times_[0] - 1.0e-14 && t_equivalent <= times_[times_.size() - 1] + 1.0e-14)
  {
    return true;
  }
  else
  {
    return false;
  }
}


CORE::UTILS::MultiFunctionVariable::MultiFunctionVariable(std::string name,
    std::vector<double> times, std::vector<std::string> description_vec,
    struct Periodicstruct periodicdata)
    : FunctionVariable(std::move(name)),
      times_(std::move(times)),
      periodic_(periodicdata.periodic),
      t1_(periodicdata.t1),
      t2_(periodicdata.t2)
{
  // create vectors of timefunction and timederivative
  timefunction_.resize(times_.size() - 1);
  for (unsigned int n = 0; n < times_.size() - 1; ++n)
  {
    timefunction_[n] =
        Teuchos::rcp(new CORE::UTILS::SymbolicExpression<double>(description_vec[n]));
  }
}


double CORE::UTILS::MultiFunctionVariable::Value(const double t)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t >= t1_ - 1.0e-14 and t <= t2_ + 1.0e-14)
  {
    t_equivalent = fmod(t + 1.0e-14, times_[times_.size() - 1] - times_[0]) - 1.0e-14;
  }

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t_equivalent < times_[0] - 1.0e-14)
  {
    dserror("t_equivalent is smaller than the starting time");
    return 0.0;
  }
  else
  {
    while (t_temp < t_equivalent - 1.0e-14)
    {
      index++;
      t_temp = times_[index];
      if (index == times_.size())
      {
        dserror("t_equivalent is bigger than the ending time");
        return 0.0;
      }
    }
  }

  // evaluate the value of the function considering the different possibilities
  double value;
  if (index == 0)
  {
    value = timefunction_[0]->Value({{"t", t_equivalent}});
  }
  else
  {
    value = timefunction_[index - 1]->Value({{"t", t_equivalent}});
  }

  return value;
}


double CORE::UTILS::MultiFunctionVariable::TimeDerivativeValue(const double t, const unsigned deg)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t >= t1_ - 1.0e-14 and t <= t2_ + 1.0e-14)
  {
    t_equivalent = fmod(t + 1.0e-14, times_[times_.size() - 1] - times_[0]) - 1.0e-14;
  }

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t_equivalent < times_[0])
  {
    dserror("t_equivalent is smaller than the starting time");
    return 0.0;
  }
  else
  {
    while (t_temp < t_equivalent - 1.0e-14)
    {
      index++;
      t_temp = times_[index];
      if (index == times_.size())
      {
        dserror("t_equivalent is bigger than the ending time");
        return 0.0;
      }
    }
  }

  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> tfad(1, 0, t_equivalent);
  tfad.val() = Sacado::Fad::DFad<double>(1, 0, t_equivalent);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> vfad;

  // evaluate the derivative of the function considering the different possibilities
  if (index == 0)
  {
    vfad = timefunction_[0]->SecondDerivative({{"t", tfad}}, {});
  }
  else
  {
    vfad = timefunction_[index - 1]->SecondDerivative({{"t", tfad}}, {});
  }

  if (deg == 0)
  {
    const double v = Value(t);
    return v;
  }
  if (deg == 1)
  {
    const double dv_dt = vfad.dx(0).val();
    return dv_dt;
  }
  else if (deg == 2)
  {
    const double d2v_dt2 = vfad.dx(0).dx(0);
    return d2v_dt2;
  }
  else
  {
    dserror("Higher than second derivative is not implemented!");
    return 0.0;
  }
}


bool CORE::UTILS::MultiFunctionVariable::ContainTime(const double t)
{
  /// check the inclusion of the considered time
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t >= t1_ - 1.0e-14 and t <= t2_ + 1.0e-14)
  {
    t_equivalent = fmod(t + 1.0e-14, times_[times_.size() - 1] - times_[0]) - 1.0e-14;
  }

  if (t_equivalent >= times_[0] - 1.0e-14 && t_equivalent <= times_[times_.size() - 1] + 1.0e-14)
  {
    return true;
  }
  else
  {
    return false;
  }
}


CORE::UTILS::FourierInterpolationVariable::FourierInterpolationVariable(std::string name,
    std::vector<double> times, std::vector<double> values, struct Periodicstruct periodicdata)
    : FunctionVariable(std::move(name)),
      times_(std::move(times)),
      values_(std::move(values)),
      periodic_(periodicdata.periodic),
      t1_(periodicdata.t1),
      t2_(periodicdata.t2)
{
}

double CORE::UTILS::FourierInterpolationVariable::Value(const double t)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t >= t1_ - 1.0e-14 and t <= t2_ + 1.0e-14)
  {
    t_equivalent = fmod(t + 1.0e-14, times_[times_.size() - 1] - times_[0]) - 1.0e-14;
  }

  const auto value = Value<double>(t_equivalent);

  return value;
}

template <typename ScalarT>
ScalarT CORE::UTILS::FourierInterpolationVariable::Value(const ScalarT& t)
{
  // source: https://en.wikipedia.org/wiki/Trigonometric_interpolation
  ScalarT value = 0.0;

  // number of interpolation nodes
  auto N = (const double)times_.size();

  // adjusting the spacing of the given independent variable
  const double scale = (times_[1] - times_[0]) * N / 2;

  // evaluate interpolant
  for (unsigned int k = 1; k <= times_.size(); ++k)
  {
    ScalarT tau;
    ScalarT xt = (t - times_[k - 1]) / scale;
    if (xt >= -1.0e-14 and xt <= 1.0e-14)
    {
      tau = 1;
    }
    else
    {
      const int mod = (int)N % 2;
      if (mod == 1)  // odd
      {
        tau = sin(N * M_PI * xt / 2) / (N * sin(M_PI * xt / 2));
      }
      else  // even
      {
        tau = sin(N * M_PI * xt / 2) / (N * tan(M_PI * xt / 2));
      }
    }

    value += values_[k - 1] * tau;
  }

  return value;
}


double CORE::UTILS::FourierInterpolationVariable::TimeDerivativeValue(
    const double t, const unsigned deg)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t >= t1_ - 1.0e-14 and t <= t2_ + 1.0e-14)
  {
    t_equivalent = fmod(t + 1.0e-14, times_[times_.size() - 1] - times_[0]) - 1.0e-14;
  }

  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> tfad(1, 0, t_equivalent);
  tfad.val() = Sacado::Fad::DFad<double>(1, 0, t_equivalent);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> vfad =
      Value<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>(tfad);

  if (deg == 0)
  {
    const double v = Value(t);
    return v;
  }
  if (deg == 1)
  {
    const double dv_dt = vfad.dx(0).val();
    return dv_dt;
  }
  else if (deg == 2)
  {
    const double d2v_dt2 = vfad.dx(0).dx(0);
    return d2v_dt2;
  }
  else
  {
    dserror("Higher than second derivative is not implemented!");
    return 0.0;
  }
}


bool CORE::UTILS::FourierInterpolationVariable::ContainTime(const double t)
{
  /// check the inclusion of the considered time
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t >= t1_ - 1.0e-14 and t <= t2_ + 1.0e-14)
  {
    t_equivalent = fmod(t + 1.0e-14, times_[times_.size() - 1] - times_[0]) - 1.0e-14;
  }

  if (t_equivalent >= times_[0] - 1.0e-14 && t_equivalent <= times_[times_.size() - 1] + 1.0e-14)
  {
    return true;
  }
  else
  {
    return false;
  }
}


CORE::UTILS::PiecewiseVariable::PiecewiseVariable(
    const std::string& name, std::vector<Teuchos::RCP<FunctionVariable>> pieces)
    : FunctionVariable(name), pieces_(std::move(pieces))
{
  if (pieces_.empty())
    dserror("A PiecewiseVariable must have at least one FunctionVariable piece.");
}


double CORE::UTILS::PiecewiseVariable::Value(const double t)
{
  return FindPieceForTime(t).Value(t);
}


double CORE::UTILS::PiecewiseVariable::TimeDerivativeValue(const double t, const unsigned int deg)
{
  return FindPieceForTime(t).TimeDerivativeValue(t, deg);
}


bool CORE::UTILS::PiecewiseVariable::ContainTime(const double t)
{
  const auto active_piece =
      std::find_if(pieces_.begin(), pieces_.end(), [t](auto& var) { return var->ContainTime(t); });

  return active_piece != pieces_.end();
}


CORE::UTILS::FunctionVariable& CORE::UTILS::PiecewiseVariable::FindPieceForTime(const double t)
{
  auto active_piece =
      std::find_if(pieces_.begin(), pieces_.end(), [t](auto& var) { return var->ContainTime(t); });

  if (active_piece == pieces_.end())
    dserror("Piece-wise variable <%s> is not defined at time %f.", Name().c_str(), t);

  return **active_piece;
}

std::vector<double> CORE::UTILS::INTERNAL::ExtractTimeVector(const INPUT::LineDefinition& timevar)
{
  // read the number of points
  int numpoints;
  timevar.ExtractInt("NUMPOINTS", numpoints);

  // read whether times are defined by number of points or by vector
  bool bynum = timevar.HasString("BYNUM");

  // read respectively create times vector
  std::vector<double> times = std::invoke(
      [&]()
      {
        if (bynum)  // times defined by number of points
        {
          // read the time range
          std::vector<double> timerange;
          timevar.ExtractDoubleVector("TIMERANGE", timerange);

          std::vector<double> times;

          // get initial and final time
          double t_initial = timerange[0];
          double t_final = timerange[1];

          // build the vector of times
          times.push_back(t_initial);
          int n = 0;
          double dt = (t_final - t_initial) / (numpoints - 1);
          while (times[n] + dt <= t_final + 1.0e-14)
          {
            if (times[n] + 2 * dt <= t_final + 1.0e-14)
            {
              times.push_back(times[n] + dt);
            }
            else
            {
              times.push_back(t_final);
            }
            ++n;
          }
          return times;
        }
        else  // times defined by vector
        {
          std::vector<double> times;
          timevar.ExtractDoubleVector("TIMES", times);
          return times;
        }
      });

  // check if the times are in ascending order
  if (!std::is_sorted(times.begin(), times.end())) dserror("the TIMES must be in ascending order");

  return times;
}

FOUR_C_NAMESPACE_CLOSE
