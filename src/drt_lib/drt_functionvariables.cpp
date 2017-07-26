/*----------------------------------------------------------------------*/
/*!
\file drt_functionvariables.cpp

\brief Time dependent variables for function definition

<pre>
\level 0

\maintainer Andrea La Spina
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#include <Sacado.hpp>
#include "drt_functionvariables.H"
#include "drt_parser.H"


DRT::UTILS::FunctionVariable::FunctionVariable(std::string name):
  name_(name)
{
}


DRT::UTILS::ParsedFunctionVariable::ParsedFunctionVariable(std::string name, std::string buf):
    FunctionVariable(name),
    timefunction_(Teuchos::rcp(new DRT::PARSER::Parser<double>(buf))),
    timederivative_(Teuchos::rcp(new DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > >(buf)))
{
  // set the timefunction
  timefunction_->AddVariable("t",0);
  timefunction_->ParseFunction();

  // set the timederivative
  timederivative_->AddVariable("t",0);
  timederivative_->ParseFunction();
}


double DRT::UTILS::ParsedFunctionVariable::Value(const double t)
{
  // set the time variable
  timefunction_->SetValue("t",t);

  // evaluate the value of the function
  double value = timefunction_->Evaluate();

  return value;
}


double DRT::UTILS::ParsedFunctionVariable::TimeDerivativeValue(const double t, const unsigned deg)
{
  Sacado::Fad::DFad< Sacado::Fad::DFad<double> > tfad(1, 0, t);
  tfad.val() = Sacado::Fad::DFad<double>(1, 0, t);
  Sacado::Fad::DFad< Sacado::Fad::DFad<double> > vfad;

  timederivative_->SetValue("t",tfad);
  vfad = timederivative_->Evaluate();

  if (deg==0)
  {
    const double v = Value(t);
    return v;
  }
  if (deg==1)
  {
    const double dv_dt = vfad.dx(0).val();
    return dv_dt;
  }
  else if (deg==2)
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


bool DRT::UTILS::ParsedFunctionVariable::ContainTime(const double t)
{
  return true;
}


DRT::UTILS::LinearInterpolationVariable::LinearInterpolationVariable(std::string name, std::vector<double> times, std::vector<double> values, struct periodicstruct periodicdata):
    FunctionVariable(name),
    times_(times),
    values_(values),
    periodic_(periodicdata.periodic),
    t1_(periodicdata.t1),
    t2_(periodicdata.t2)
{
}


double DRT::UTILS::LinearInterpolationVariable::Value(const double t)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t>=t1_-1.0e-14 and t<=t2_+1.0e-14)
  {
    t_equivalent = fmod(t+1.0e-14,times_[times_.size()-1] - times_[0])-1.0e-14;
  }

  const double value = Value<double>(t_equivalent);

  return value;
}

template <typename ScalarT>
ScalarT DRT::UTILS::LinearInterpolationVariable::Value(const ScalarT& t)
{
  ScalarT value = 0.0;

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t<times_[0]-1.0e-14)
  {
    dserror("t_equivalent is smaller than the starting time");
    return 0.0;
  }
  else if (t<=times_[0]+1.0e-14)
  {
    index=1;
  }
  else
  {
    while(t_temp<=t-1.0e-14)
    {
      index++;
      t_temp=times_[index];

      if(index==times_.size())
      {
        dserror("t_equivalent is bigger than the ending time");
        return 0.0;
      }
    }
  }

  // evaluate the value of the function
  value = values_[index-1] + (values_[index]-values_[index-1])/(times_[index]-times_[index-1]) * (t-times_[index-1]);

  return value;
}


double DRT::UTILS::LinearInterpolationVariable::TimeDerivativeValue(const double t, const unsigned deg)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t>=t1_-1.0e-14 and t<=t2_+1.0e-14)
  {
    t_equivalent = fmod(t+1.0e-14,times_[times_.size()-1] - times_[0])-1.0e-14;
  }

  Sacado::Fad::DFad< Sacado::Fad::DFad<double> > tfad(1, 0, t_equivalent);
  tfad.val() = Sacado::Fad::DFad<double>(1, 0, t_equivalent);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double> > vfad = Value< Sacado::Fad::DFad< Sacado::Fad::DFad<double> > >(tfad);

  if (deg==0)
  {
    const double v = Value(t);
    return v;
  }
  if (deg==1)
  {
    const double dv_dt = vfad.dx(0).val();
    return dv_dt;
  }
  else if (deg==2)
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


bool DRT::UTILS::LinearInterpolationVariable::ContainTime(const double t)
{
  /// check the inclusion of the considered time
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t>=t1_-1.0e-14 and t<=t2_+1.0e-14)
  {
    t_equivalent = fmod(t+1.0e-14,times_[times_.size()-1] - times_[0])-1.0e-14;
  }

  if (t_equivalent >= times_[0]-1.0e-14 && t_equivalent <= times_[times_.size()-1]+1.0e-14)
  {
    return true;
  }
  else
  {
    return false;
  }
}


DRT::UTILS::MultiFunctionVariable::MultiFunctionVariable(std::string name, std::vector<double> times, std::vector<std::string> description_vec, struct periodicstruct periodicdata):
    FunctionVariable(name),
    times_(times),
    periodic_(periodicdata.periodic),
    t1_(periodicdata.t1),
    t2_(periodicdata.t2)
{
  // create vectors of timefunction and timederivative
  timefunction_.resize(times_.size()-1);
  timederivative_.resize(times_.size()-1);
  for (unsigned int n=0; n<times_.size()-1; ++n)
  {
    timefunction_[n] = Teuchos::rcp(new DRT::PARSER::Parser<double>(description_vec[n]));
    timefunction_[n]->AddVariable("t",0);
    timefunction_[n]->ParseFunction();

    timederivative_[n] = Teuchos::rcp(new DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > >(description_vec[n]));
    timederivative_[n]->AddVariable("t",0);
    timederivative_[n]->ParseFunction();
  }
}


double DRT::UTILS::MultiFunctionVariable::Value(const double t)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t>=t1_-1.0e-14 and t<=t2_+1.0e-14)
  {
    t_equivalent = fmod(t+1.0e-14,times_[times_.size()-1] - times_[0])-1.0e-14;
  }

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t_equivalent<times_[0]-1.0e-14)
  {
    dserror("t_equivalent is smaller than the starting time");
    return 0.0;
  }
  else
  {
    while(t_temp<t_equivalent-1.0e-14)
    {
      index++;
      t_temp=times_[index];
      if(index==times_.size())
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
    timefunction_[0]->SetValue("t",t_equivalent);
    value = timefunction_[0]->Evaluate();
  }
  else
  {
    timefunction_[index-1]->SetValue("t",t_equivalent);
    value = timefunction_[index-1]->Evaluate();
  }

  return value;
}


double DRT::UTILS::MultiFunctionVariable::TimeDerivativeValue(const double t, const unsigned deg)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t>=t1_-1.0e-14 and t<=t2_+1.0e-14)
  {
    t_equivalent = fmod(t+1.0e-14,times_[times_.size()-1] - times_[0])-1.0e-14;
  }

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t_equivalent<times_[0])
  {
    dserror("t_equivalent is smaller than the starting time");
    return 0.0;
  }
  else
  {
    while(t_temp<t_equivalent-1.0e-14)
    {
      index++;
      t_temp=times_[index];
      if(index==times_.size())
      {
        dserror("t_equivalent is bigger than the ending time");
        return 0.0;
      }
    }
  }

  Sacado::Fad::DFad< Sacado::Fad::DFad<double> > tfad(1, 0, t_equivalent);
  tfad.val() = Sacado::Fad::DFad<double>(1, 0, t_equivalent);
  Sacado::Fad::DFad< Sacado::Fad::DFad<double> > vfad;

  // evaluate the derivative of the function considering the different possibilities
  if (index == 0)
  {
    timederivative_[0]->SetValue("t",tfad);
    vfad = timederivative_[0]->Evaluate();
  }
  else
  {
    timederivative_[index-1]->SetValue("t",tfad);
    vfad = timederivative_[index-1]->Evaluate();
  }

  if (deg==0)
  {
    const double v = Value(t);
    return v;
  }
  if (deg==1)
  {
    const double dv_dt = vfad.dx(0).val();
    return dv_dt;
  }
  else if (deg==2)
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


bool DRT::UTILS::MultiFunctionVariable::ContainTime(const double t)
{

  /// check the inclusion of the considered time
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t>=t1_-1.0e-14 and t<=t2_+1.0e-14)
  {
    t_equivalent = fmod(t+1.0e-14,times_[times_.size()-1] - times_[0])-1.0e-14;
  }

  if (t_equivalent >= times_[0]-1.0e-14 && t_equivalent <= times_[times_.size()-1]+1.0e-14)
  {
    return true;
  }
  else
  {
    return false;
  }
}


DRT::UTILS::FourierInterpolationVariable::FourierInterpolationVariable(std::string name, std::vector<double> times, std::vector<double> values, struct periodicstruct periodicdata):
    FunctionVariable(name),
    times_(times),
    values_(values),
    periodic_(periodicdata.periodic),
    t1_(periodicdata.t1),
    t2_(periodicdata.t2)
{
}

double DRT::UTILS::FourierInterpolationVariable::Value(const double t)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t>=t1_-1.0e-14 and t<=t2_+1.0e-14)
  {
    t_equivalent = fmod(t+1.0e-14,times_[times_.size()-1] - times_[0])-1.0e-14;
  }

  const double value = Value<double>(t_equivalent);

  return value;
}

template <typename ScalarT>
ScalarT DRT::UTILS::FourierInterpolationVariable::Value(const ScalarT& t)
{
  // source: https://en.wikipedia.org/wiki/Trigonometric_interpolation
  ScalarT value = 0.0;
  const double PI = 3.14159265358979323846;

  // number of interpolation nodes
  double N = (const double)times_.size();

  // adjusting the spacing of the given independent variable
  const double scale = (times_[1] -times_[0]) *N / 2;

  // evaluate interpolant
  for (unsigned int k=1; k<=times_.size(); ++k)
  {
    ScalarT tau;
    ScalarT xt = (t - times_[k-1])/scale;
    if (xt >= -1.0e-14 and xt <= 1.0e-14)
    {
      tau = 1;
    }
    else
    {
      const int mod = (int)N %2;
      if (mod == 1) // odd
      {
        tau = sin(N*PI*xt/2) / (N*sin(PI*xt/2));
      }
      else          // even
      {
        tau = sin(N*PI*xt/2) / (N*tan(PI*xt/2));
      }
    }

    value += values_[k-1] * tau;
  }

  return value;
}


double DRT::UTILS::FourierInterpolationVariable::TimeDerivativeValue(const double t, const unsigned deg)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t>=t1_-1.0e-14 and t<=t2_+1.0e-14)
  {
    t_equivalent = fmod(t+1.0e-14,times_[times_.size()-1] - times_[0])-1.0e-14;
  }

  Sacado::Fad::DFad< Sacado::Fad::DFad<double> > tfad(1, 0, t_equivalent);
  tfad.val() = Sacado::Fad::DFad<double>(1, 0, t_equivalent);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double> > vfad = Value< Sacado::Fad::DFad< Sacado::Fad::DFad<double> > >(tfad);

  if (deg==0)
  {
    const double v = Value(t);
    return v;
  }
  if (deg==1)
  {
    const double dv_dt = vfad.dx(0).val();
    return dv_dt;
  }
  else if (deg==2)
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


bool DRT::UTILS::FourierInterpolationVariable::ContainTime(const double t)
{
  /// check the inclusion of the considered time
  double t_equivalent = t;

  // handle periodicity
  if (periodic_ and t>=t1_-1.0e-14 and t<=t2_+1.0e-14)
  {
    t_equivalent = fmod(t+1.0e-14,times_[times_.size()-1] - times_[0])-1.0e-14;
  }

  if (t_equivalent >= times_[0]-1.0e-14 && t_equivalent <= times_[times_.size()-1]+1.0e-14)
  {
    return true;
  }
  else
  {
    return false;
  }
}
