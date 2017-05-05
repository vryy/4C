/*----------------------------------------------------------------------*/
/*!
\file drt_functionvariables.cpp

\brief Time dependent variables for function definition

<pre>
\level 1

\maintainer Andrea La Spina
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#include <Sacado.hpp>
#include "drt_functionvariables.H"
#include "drt_parser.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ParsedFunctionVariable::~ParsedFunctionVariable()
{}


DRT::UTILS::ParsedFunctionVariable::ParsedFunctionVariable(std::string name, std::string buf)
{
  // set the name
  name_ = name;

  // set the timefunction
  timefunction_ = Teuchos::rcp(new DRT::PARSER::Parser<double>(buf));
  timefunction_->AddVariable("t",0);
  timefunction_->ParseFunction();

  // set the timederivative
  timederivative_ = Teuchos::rcp(new DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > >(buf));
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

  double der;
  if (deg==1)
  {
    double dv_dt = vfad.dx(0).val();
    der = dv_dt;
    return der;
  }
  else if (deg==2)
  {
    double d2v_dt2 = vfad.dx(0).dx(0);
    der = d2v_dt2;
    return der;
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
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::LinearInterpolationVariable::~LinearInterpolationVariable()
{}


DRT::UTILS::LinearInterpolationVariable::LinearInterpolationVariable(std::string name, std::vector<double> times, std::vector<double> values, struct periodicstruct periodicdata)
{
  // set the name
  name_ = name;

  // set the times
  times_ = times;

  // set the values
  values_ = values;

  // set the periodic data
  periodic_ = periodicdata.periodic;
  t1_ = periodicdata.t1;
  t2_ = periodicdata.t2;
}


double DRT::UTILS::LinearInterpolationVariable::Value(const double t)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;
  if (periodic_)
  {
    t_equivalent = EquivalentTime(t);
  }

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t_equivalent<times_[0])
    return 0;
  else if (t_equivalent==times_[0])
  {
    return values_[0];
  }
  else
  {
    while(t_temp<t_equivalent-1.0e-13)
    {
      index++;
      t_temp=times_[index];

      if(index==times_.size())
        return 0;
    }
  }
  // evaluate the value of the function
  double value = values_[index-1] + (values_[index]-values_[index-1])/(times_[index]-times_[index-1]) * (t_equivalent-times_[index-1]);

  return value;
}


double DRT::UTILS::LinearInterpolationVariable::TimeDerivativeValue(const double t, const unsigned deg)
{
  // evaluate an equivalent time for a periodic variable
  double t_original = t;
  double t_equivalent = t;
  bool shift = false;
  if (periodic_)
  {
    t_equivalent = EquivalentTime(t);
    // check if the time has been shifted due to periodicity
    if(!(t_equivalent<t_original+1.0e-13 && t_equivalent>t_original-1.0e-13))
      shift = true;
  }

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t_equivalent<times_[0])
    return 0;
  else
  {
    while(t_temp<t_equivalent-1.0e-13)
    {
      index++;
      t_temp=times_[index];
      if(index==times_.size())
        return 0;
    }
  }

  // evaluate the derivative of the function considering the different possibilities
  double der;
  if (deg==1)
  {
    if (shift) // evaluate the derivative inside a periodic repetition
    {
      if (t_equivalent != times_[index]) // evaluate the derivative inside a slice
      {
        der = (values_[index]-values_[index-1])/(times_[index]-times_[index-1]);
      }
      else // evaluate derivative in a sampling time
      {
        if (index == 0) // if the time is the first one
        {
          if (t_original == t1_) // if the original time was the initial time of repetition
          {
            if (t1_ == times_[times_.size()-1])
            {
              // if the initial time of repetition coincides with the final time of variable
              // consider the average value of the initial and final derivatives
              double der_initial;
              double der_final;

              der_initial = (values_[1] - values_[0]) / (times_[1] - times_[0]);
              der_final   = (values_[times_.size()-1] - values_[times_.size()-2]) / (times_[times_.size()-1] - times_[times_.size()-2]);
              der = (der_initial + der_final) / 2;
            }
            else
            {
              // if the initial time of repetition does not coincide with the final time of variable
              // consider the initial derivative
              der = (values_[1] - values_[0]) / (times_[1] - times_[0]);
            }
          }
          else if (t_original == t2_) // if the original time was the final time of repetition
          {
            // consider the final derivative
            der = (values_[times_.size()-1] - values_[times_.size()-2]) / (times_[times_.size()-1] - times_[times_.size()-2]);
          }
          else // if the original time was not the initial or final time of repetition
          {
            // consider the average value of the initial and final derivatives
            double der_initial;
            double der_final;

            der_initial = (values_[1] - values_[0]) / (times_[1] - times_[0]);
            der_final   = (values_[times_.size()-1] - values_[times_.size()-2]) / (times_[times_.size()-1] - times_[times_.size()-2]);
            der = (der_initial + der_final) / 2;
          }
        }
        else if (index == times_.size() - 1) // if the time is the last one
        {
          if (t_original == t2_) // if the original time was the final time of repetition
          {
            // consider the final derivative
            der = (values_[times_.size()-1] - values_[times_.size()-2]) / (times_[times_.size()-1] - times_[times_.size()-2]);
          }
          else  // if the original time was not the final time of repetition
          {
            // consider the average value of the initial and final derivatives
            double der_initial;
            double der_final;

            der_initial = (values_[1] - values_[0]) / (times_[1] - times_[0]);
            der_final   = (values_[times_.size()-1] - values_[times_.size()-2]) / (times_[times_.size()-1] - times_[times_.size()-2]);
            der = (der_initial + der_final) / 2;
          }
        }
        else // if the time is an internal one take the average value of the derivatives
        {
          double der_left;
          double der_right;

          der_left =  (values_[index] - values_[index-1]) / (times_[index] - times_[index-1]);
          der_right = (values_[index+1] - values_[index]) / (times_[index+1] - times_[index]);
          der = (der_left + der_right) / 2;
        }
      }
    }
    else // evaluate the derivative outside a (possible) periodic repetition
    {
      if (t_equivalent != times_[index]) // evaluate the derivative inside a slice
      {
        der = (values_[index]-values_[index-1])/(times_[index]-times_[index-1]);
      }
      else // evaluate derivative in a sampling time
      {
        if (index == 0) // if the time is the first one consider the derivative on the right
        {
          der = (values_[index+1] - values_[index]) / (times_[index+1] - times_[index]);
        }
        else if (index == times_.size() - 1) // if the time is the last one consider the derivative on the left
        {
          der = (values_[index] - values_[index-1]) / (times_[index] - times_[index-1]);
        }
        else // if the time is an internal one take the average value of the derivatives
        {
          double der_left;
          double der_right;

          der_left =  (values_[index] - values_[index-1]) / (times_[index] - times_[index-1]);
          der_right = (values_[index+1] - values_[index]) / (times_[index+1] - times_[index]);
          der = (der_left + der_right) / 2;
        }
      }
    }
  }
  else
  {
    der = 0;
  }

  return der;
}


double DRT::UTILS::LinearInterpolationVariable::EquivalentTime(const double t)
{
  double t_equivalent;

  if (t>=times_[0] && t<=times_[times_.size()-1]+1E-13)
  {
    return t;
  }
  else if (t>=t1_ && t<=t2_)
  {
    // evaluate the period
    double T = times_[times_.size()-1] - times_[0];

    // evaluate the number of periods to subtract
    double n = floor ((t - t1_) / T);

    // evaluate the equivalent time
    t_equivalent = times_[0] + ((t - n * T) - t1_);

    return t_equivalent;
  }
  else
  {
    return t;
  }
}


bool DRT::UTILS::LinearInterpolationVariable::ContainTime(const double t)
{
  /// check the inclusion of the considered time
  double t_equivalent = t;
  if (periodic_)
  {
    t_equivalent = EquivalentTime(t);
  }

  if (times_[0] <= t_equivalent && t_equivalent <= times_[times_.size()-1]+1.0e-12)
  {
    return true;
  }
  else
  {
    return false;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::MultiFunctionVariable::~MultiFunctionVariable()
{}


DRT::UTILS::MultiFunctionVariable::MultiFunctionVariable(std::string name, std::vector<double> times, std::vector<std::string> description_vec, struct periodicstruct periodicdata)
{
  // set the name
  name_ = name;

  // set the times
  times_ = times;

  // set the periodic data
  periodic_ = periodicdata.periodic;
  t1_ = periodicdata.t1;
  t2_ = periodicdata.t2;

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
  if (periodic_)
  {
    t_equivalent = EquivalentTime(t);
  }

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t_equivalent<times_[0])
    return 0;
  else
  {
    while(t_temp<t_equivalent-1.0e-13)
    {
      index++;
      t_temp=times_[index];
      if(index==times_.size())
        return 0;
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
  if (periodic_)
  {
    t_equivalent = EquivalentTime(t);
  }

  // find the right time slice
  double t_temp = times_[0];
  unsigned int index = 0;
  if (t_equivalent<times_[0])
    return 0;
  else
  {
    while(t_temp<t_equivalent-1.0e-13)
    {
      index++;
      t_temp=times_[index];
      if(index==times_.size())
        return 0;
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

  double der;
  if (deg==1)
  {
    double dv_dt = vfad.dx(0).val();
    der = dv_dt;
    return der;
  }
  else if (deg==2)
  {
    double d2v_dt2 = vfad.dx(0).dx(0);
    der = d2v_dt2;
    return der;
  }
  else
  {
    dserror("Higher than second derivative is not implemented!");
    return 0.0;
  }
}


double DRT::UTILS::MultiFunctionVariable::EquivalentTime(const double t)
{
  double t_equivalent;

  if (t>=times_[0] && t<=times_[times_.size()-1]+1E-13)
  {
    return t;
  }
  else if (t>=t1_ && t<=t2_)
  {
    // evaluate the period
    double T = times_[times_.size()-1] - times_[0];

    // evaluate the number of periods to subtract
    double n = floor ((t - t1_) / T);

    // evaluate the equivalent time
    t_equivalent = times_[0] + ((t - n * T) - t1_);

    return t_equivalent;
  }
  else
  {
    return t;
  }
}


bool DRT::UTILS::MultiFunctionVariable::ContainTime(const double t)
{

  /// check the inclusion of the considered time
  double t_equivalent = t;
  if (periodic_)
  {
    t_equivalent = EquivalentTime(t);
  }

  if (times_[0] <= t_equivalent && t_equivalent <= times_[times_.size()-1]+1.0e-13)
  {
    return true;
  }
  else
  {
    return false;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::FourierInterpolationVariable::~FourierInterpolationVariable()
{}


DRT::UTILS::FourierInterpolationVariable::FourierInterpolationVariable(std::string name, std::vector<double> times, std::vector<double> values, struct periodicstruct periodicdata)
{
  // set the name
  name_ = name;

  // set the times
  times_ = times;

  // set the values
  values_ = values;

  // set the periodic data
  periodic_ = periodicdata.periodic;
  t1_ = periodicdata.t1;
  t2_ = periodicdata.t2;
}

double DRT::UTILS::FourierInterpolationVariable::Value(const double t)
{
  const double value = Value<double>(t);

  return value;
}

template <typename ScalarT>
ScalarT DRT::UTILS::FourierInterpolationVariable::Value(const ScalarT& t)
{
  const int DataLength=values_.size();
  std::vector<double> SampleNumber;
  std::vector<double> EvenCoefficient;
  std::vector<double> OddCoefficient;
  ScalarT fac = 0.0;
  double C = (const double)values_.size();
  double PI = 3.14159265358979323846;
  double period = times_[times_.size()-1] - times_[0];

  for (int p=0; p<DataLength; p++)
  {
    SampleNumber.push_back(values_[p]);
  }

  for (int p=0; p<=DataLength/2; p++)
  {
    EvenCoefficient.push_back(0);
    OddCoefficient.push_back(0);
  }

  for (int p=0; p<=DataLength/2; p++)
  {
    EvenCoefficient[p] = 0;
    OddCoefficient[p] = 0;

    for (int num=0; num<=DataLength-1; num++)
    {
      EvenCoefficient[p] = EvenCoefficient[p]
                                           + 2/C*SampleNumber[num]*cos(2*PI*p*(num+1)/C);
      OddCoefficient[p] = OddCoefficient[p]
                                         + 2/C*SampleNumber[num]*sin(2*PI*p*(num+1)/C);
    }
  }

  EvenCoefficient[DataLength/2] = EvenCoefficient[DataLength/2]/2;
  OddCoefficient[DataLength/2] = 0;
  fac = EvenCoefficient[0]/2;

  for (int h=1; h<=DataLength/2; h++)
  {
    fac = fac
        + (ScalarT)EvenCoefficient[h]*cos(2*PI*h*t/period)
        + (ScalarT)OddCoefficient[h]*sin(2*PI*h*t/period);
  }

  return fac;
}

// The implementation of the Fourier interpolation produces a wrong shift of the function to the right.
// Here there is the code I have implemented for a correct use of the Fourier interpolation.
// Andrea La Spina 01/02/2017

//double DRT::UTILS::FourierInterpolationVariable::Value(const double t)
//{
//  // MY CODE **************************************************************************************
//  // evaluate an equivalent time for a periodic variable
//  double t_equivalent = t;
//  if (periodic_)
//  {
//    t_equivalent = EquivalentTime(t);
//  }
//
//  // source: https://en.wikipedia.org/wiki/Trigonometric_interpolation
//  double PI = 3.14159265358979323846;
//
//  // number of interpolation nodes
//  double N = (const double)times_.size();
//
//  // adjusting the spacing of the given independent variable
//  double h = 2/N;
//  double scale = (times_[1] -times_[0]) / h;
//  std::vector<double> x;
//  for (unsigned int n=0; n<times_.size(); ++n)
//  {
//    x.push_back(times_[n]/scale);
//  }
//  double xi = t_equivalent / scale;
//
//  // evaluate interpolant
//  double value = 0;
//  for (int k=1; k<=N; ++k)
//  {
//    double tau;
//    double xt = xi - x[k-1];
//    if (xt == 0)
//    {
//      tau = 1;
//    }
//    else
//    {
//      int mod = (int)N %2;
//      if (mod == 1) // odd
//      {
//        tau = sin(N*PI*xt/2) / (N*sin(PI*xt/2));
//      }
//      else          // even
//      {
//        tau = sin(N*PI*xt/2) / (N*tan(PI*xt/2));
//      }
//    }
//
//    value = value + values_[k-1] * tau;
//  }
//
//  return value;
//  // **********************************************************************************************
//}


double DRT::UTILS::FourierInterpolationVariable::TimeDerivativeValue(const double t, const unsigned deg)
{
  // evaluate an equivalent time for a periodic variable
  double t_equivalent = t;
  if (periodic_)
  {
    t_equivalent = EquivalentTime(t);
  }

  Sacado::Fad::DFad< Sacado::Fad::DFad<double> > tfad(1, 0, t_equivalent);
  tfad.val() = Sacado::Fad::DFad<double>(1, 0, t_equivalent);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double> > vfad = Value< Sacado::Fad::DFad< Sacado::Fad::DFad<double> > >(tfad);

  double der;
  if (deg==1)
  {
    double dv_dt = vfad.dx(0).val();
    der = dv_dt;
    return der;
  }
  else if (deg==2)
  {
    double d2v_dt2 = vfad.dx(0).dx(0);
    der = d2v_dt2;
    return der;
  }
  else
  {
    dserror("Higher than second derivative is not implemented!");
    return 0.0;
  }
}


double DRT::UTILS::FourierInterpolationVariable::EquivalentTime(const double t)
{
  double t_equivalent;

  if (t>=times_[0] && t<=times_[times_.size()-1]+1E-13)
  {
    return t;
  }
  else if (t>=t1_ && t<=t2_)
  {
    // evaluate the period
    double T = times_[times_.size()-1] - times_[0];

    // evaluate the number of periods to subtract
    double n = floor ((t - t1_) / T);

    // evaluate the equivalent time
    t_equivalent = times_[0] + ((t - n * T) - t1_);

    return t_equivalent;
  }
  else
  {
    return t;
  }
}


bool DRT::UTILS::FourierInterpolationVariable::ContainTime(const double t)
{
  /// check the inclusion of the considered time
  double t_equivalent = t;
  if (periodic_)
  {
    t_equivalent = EquivalentTime(t);
  }

  if (times_[0] <= t_equivalent && t_equivalent <= times_[times_.size()-1])
  {
    return true;
  }
  else
  {
    return false;
  }
}
