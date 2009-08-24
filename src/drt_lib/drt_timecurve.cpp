/*----------------------------------------------------------------------*/
/*!
\file drt_timecurve.cpp

\brief Managing and evaluating of time curves

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <Sacado.hpp>

#include "drt_parser.H"
#include "drt_timecurve.H"
#include "drt_dserror.H"

#include "drt_globalproblem.H"
#include "../drt_mat/newtonianfluid.H"

#ifdef PARALLEL
#include <mpi.h>
#endif


using namespace std;
using namespace Teuchos;


namespace DRT {
namespace UTILS {

  /*====================================================================*/
  /// base of all slices of a time curve with start and end time
  class TimeSlice
  {
  public:

    /// construct with begin and end time
    TimeSlice(double begin, double end) : begin_(begin), end_(end) {}

    /// virtual destructor because this class will be derived
    virtual ~TimeSlice() {}

    /// interval check [\c begin_ , \c end_ [
    bool contains(double t) const { return t>=begin_-1e-13 and t<end_+1e-13; }

    /// evaluate time curve \c f at given time \c t
    virtual double f(double time) = 0;

    /// \brief Evaluate time curve and its derivatives
    ///
    /// Delivers vector containing at position i of resulting vector
    /// the i-th derivative of the time curve function at time \c t.
    /// The highest degree of resulting
    /// time derivatives is set by the argument \c deg.
    /// \author bborn
    /// \date 02/08
    virtual std::vector<double> FctDer(
      const double t,  ///< evaluation time
      const unsigned deg  ///< highest degree of differentiation
    ) = 0;

    /// begin time of slice
    double begin() const { return begin_; }

    /// end time of slice
    double end() const { return end_; }

    /// debug output of this slice
    virtual void Print(std::ostream& out) const = 0;

  private:
    /// begin time point of time slice
    double begin_;
    /// end time point of time slice
    double end_;

    /// output operator
    friend std::ostream& operator<<(std::ostream& out, const TimeSlice& slice);
  };


  /*--------------------------------------------------------------------*/
  /// time slice defined by polygonal (with jumps!)
  class PolygonalTimeSlice : public TimeSlice
  {
  public:
    /// constructor
    PolygonalTimeSlice(double begin, double end, double vbegin, double vend);

    /// evaluate time curve \c f at given time \c t
    double f(double t);

    /// evaluate time curve and its derivatives
    std::vector<double> FctDer(const double t, const unsigned deg);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;
  private:

    /// value at begin of time slice
    double value_begin_;

    /// value at end of time slice
    double value_end_;
  };


  /*--------------------------------------------------------------------*/
  /// explicit strings with fixed implementation
  /*!
    Depreciated. Do not extend this class. Create a new one if you
    need. Or even better use ExprTimeSlice.
  */
  class ExplicitTimeSlice : public TimeSlice
  {
  public:
    ExplicitTimeSlice(int numex, double c1, double c2);

    /// function template
    template <typename ScalarT> ScalarT Fct(const ScalarT& t);

    /// function template explicit for type double
    double Fct(const double& t);

    /// function template explicit for type Sacado::Fad::DFad< Sacado::Fad::DFad<double> >
    /// needed for automatic differentiation to get 1st & 2nd derivative
    Sacado::Fad::DFad< Sacado::Fad::DFad<double> > Fct
    (
      const Sacado::Fad::DFad< Sacado::Fad::DFad<double> >& t
    );

    /// evaluate time curve at given time
    double f(double t);

    /// evaluate time curve and its derivatives
    std::vector<double> FctDer(const double t, const unsigned deg);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;
  private:
    int numex_;
    double c1_;
    double c2_;
  };


  /*--------------------------------------------------------------------*/
  /// special respiration time slice
  class LungTimeSlice : public TimeSlice
  {
  public:
    LungTimeSlice(double frequ, double ppeep, double phase);

    /// function template
    template <typename ScalarT> ScalarT Fct(const ScalarT& t);

    /// function template explicit for type double
    double Fct(const double& t);

    /// function template explicit for type Sacado::Fad::DFad< Sacado::Fad::DFad<double> >
    /// needed for automatic differentiation to get 1st & 2nd derivative
    Sacado::Fad::DFad< Sacado::Fad::DFad<double> > Fct
    (
      const Sacado::Fad::DFad< Sacado::Fad::DFad<double> >& t
    );

    /// evaluate time curve at given time
    double f(double t);

    /// evaluate time curve and its derivatives
    std::vector<double> FctDer(const double t, const unsigned deg);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;
  private:
    /// frequency
    double frequ_;
    /// positive and expiratory pressure
    double ppeep_;
    /// phase
    double phase_;
  };

  /*--------------------------------------------------------------------*/
  /// special respiration time slice
  class BloodTimeSlice : public TimeSlice
  {
  public:
    BloodTimeSlice(double period,
                   double flowrate,
                   int points,
                   std::vector<double>& ArrayLength );

    /// function template
    template <typename ScalarT> ScalarT Fct(const ScalarT& t);

    /// function template explicit for type double
    double Fct(const double& t);

    /// function template explicit for type Sacado::Fad::DFad< Sacado::Fad::DFad<double> >
    /// needed for automatic differentiation to get 1st & 2nd derivative
    Sacado::Fad::DFad< Sacado::Fad::DFad<double> > Fct
    (
      const Sacado::Fad::DFad< Sacado::Fad::DFad<double> >& t
    );

    /// evaluate time curve at given time
    double f(double t);

    /// evaluate time curve and its derivatives
    std::vector<double> FctDer(const double t, const unsigned deg);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;
  private:
    double period_;
    double flowrate_;
    int points_;
    std::vector<double> ArrayLength_;
  };


  /*--------------------------------------------------------------------*/
  /// time slice based on parsed expression with variable t
  /*!
    A time slice that evaluates an expression at every point in
    time. The expression can be arbitrary, it is parsed by a top down
    parser. The syntax of such a curve in the input file looks like
    that:
    <pre>
    ------------------------------------------------------------CURVE1
    CURVE1 on EXPR FUNC sin(t*pi/2) t1 0.0 t2 1.0
    </pre>
    Here a one-slice curve is defined that evaluates sin(t*pi/2)
    between t=0 and t=1. The time curve would return sin(pi/2) for all
    t>1 and (theoretically) sin(0) for all t<0. But this is a feature
    of TimeCurve::f. All classes derived from TimeSlice can only be
    evaluated in the defined range. This way it is possible to define
    curves from many slices:
    <pre>
    ------------------------------------------------------------CURVE1
    CURVE1 on EXPR FUNC sin(t*pi/2) t1 0.0 t2 1.0
    CURVE1 on EXPR FUNC 3*t t1 1.0 t2 2.0
    CURVE1 on EXPR FUNC exp(2*t^(1/4)) t1 2.0 t2 3.0
    CURVE1 on EXPR FUNC acos(t/2) t1 3.0 t2 4.0
    </pre>
    But of course nobody needs such a mess.
   */
  class ExprTimeSlice : public TimeSlice
  {
  public:
    /// construct syntax tree from string buffer
    ExprTimeSlice(double begin, double end, char* buf);

    /// explicit destructor that frees the syntax tree
    ~ExprTimeSlice();

    /// evaluate time curve at given time
    double f(double t);

    /// evaluate time curve and its derivatives
    std::vector<double> FctDer(const double t, const unsigned deg);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;

  private:

    /// parsed expression as syntax tree
    DRT::PARSER::Parser<double> parsexpr_;

    /// parsed expression twice automatically differentiated
    DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > > parsexprdd_;
  };

}
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::TimeCurveManager::ReadInput()
{
  curves_.clear();

  // test for as many curves as there are
  for (int i = 1;; ++i)
  {
    ostringstream curve;
    curve << "--CURVE" << i;
    if (frfind(curve.str().c_str())==1)
    {
      frread();
      int ierr;

      // stop if there is no content to this curve
      // no further curves are read
      frchk("---",&ierr);
      if (ierr==1)
      {
        break;
      }

      int id;
      frint("CURVE",&id,&ierr);
      if (ierr!=1) dserror("cannot read CURVE%d", i);
      if (id!=i) dserror("expected CURVE%d but got CURVE%d", i, id);

      // so we have a new time curve
      curves_.push_back(TimeCurve());
      TimeCurve& curve = curves_.back();

      for (;; frread())
      {
        // we are done once we hit the border
        frchk("---",&ierr);
        if (ierr==1)
        {
          break;
        }

        frchk("Polygonal",&ierr);
        if (ierr==1)
        {
          char* colpointer = strstr(fractplace(),"BYABSTIME");
          colpointer += 13;
          double begin = strtod(colpointer,&colpointer);
          double end = strtod(colpointer,&colpointer);
          colpointer = strstr(fractplace(),"FACTOR");
          colpointer += 6;
          double vbegin = strtod(colpointer,&colpointer);
          double vend = strtod(colpointer,&colpointer);
          curve.AddSlice(rcp(new PolygonalTimeSlice(begin,end,vbegin,vend)));
          continue;
        }
        frchk("Explicit",&ierr);
        if (ierr==1)
        {
          char buffer[100];
          frchar("FUNC",buffer,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);

          int numex = 0;
          if (string(buffer)=="f(t)=sin(t:C1*PI:2)_for_t<_C1_else_f(t)=1")
            numex=-1;
          else if (string(buffer)=="f(t)=exp(1-1:t)_for_t<C1_else_const.")
            numex=-2;
          else if (string(buffer)=="f(t)=1-cos(2*PI*C1*t)")
            numex=-3;
          else if (string(buffer)=="f(t)=C2*sin(2PI*C1*t)")
            numex=-4;
          else if (string(buffer)=="f(t)=(sin(PI(t:C1-0.5))+1)*0.5")
            numex=-5;
          else if (string(buffer)=="BELTRAMI")
            numex=-6;
          else if (string(buffer)=="KIM-MOIN")
            numex=-7;
          else if (string(buffer)=="f(t)=(C2:2PI*C1)*cos(2PI*C1*t)")
            numex=-8;
          else if (string(buffer)=="f(t)=t:2-C1:(2PI)*cos(PI*t:C1-PI:2)")
            numex=-9; /* time integral of numex -5 */
          else if (string(buffer)=="f(t)=-0.5*cos(PI*(T-C1)/C2)+0.5")
            numex=-11; /* ramp function with horizontal slopes */
          else if (string(buffer)=="f(t)=(1-C2/C1)*T+C2")
            numex=-12; /* linearly increasing function with non-zero initial value */
          else
            dserror("Cannot read function of CURVE%d: %s",i,string(buffer).c_str());

          double c1;
          double c2;

          frdouble("c1",&c1,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);
          frdouble("c2",&c2,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);

          curve.AddSlice(rcp(new ExplicitTimeSlice(numex,c1,c2)));
          continue;
        }
        frchk("EXPR",&ierr);
        if (ierr==1)
        {
          char buffer[500];
          frchar("FUNC",buffer,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);

          double begin;
          double end;

          frdouble("t1",&begin,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);
          frdouble("t2",&end,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);

          curve.AddSlice(rcp(new ExprTimeSlice(begin, end, buffer)));
          continue;
        }
        frchk("LungSinus",&ierr);
        if (ierr==1)
        {
          double frequ;
          double ppeep;
          double phase;

          frdouble("Frequ",&frequ,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);
          frdouble("pPEEP",&ppeep,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);
          frdouble("Phase",&phase,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);

          curve.AddSlice(rcp(new LungTimeSlice(frequ, ppeep, phase)));
          continue;
        }
        frchk("PhysiologicalWaveform",&ierr);
        if (ierr==1)
        {
          double period;
          double flowrate;
          int points;
          std::vector<double> ArrayLength(60);

          frdouble("Period",&period,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);
          frdouble("Flowrate",&flowrate,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);
          frint("Samplingpoints",&points,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);
          frdouble_n("Arrayread",&(ArrayLength[0]),points,&ierr);
          if (ierr!=1) dserror("cannot read CURVE%d",i);


          curve.AddSlice(rcp(new BloodTimeSlice(period, flowrate, points, ArrayLength)));
          continue;
        }
        dserror("unknown type of time curve in CURVE%d", i);
      }
    }
    else
    {
      // there is no such curve, stop reading
      break;
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::TimeCurve& DRT::UTILS::TimeCurveManager::Curve(int num)
{
  // ensure that desired curve is available (prevents segmentation fault)
  if ((curves_.size()< (unsigned int)(num+1)) || num<0)
    dserror("time curve %d not available",num+1);

  return curves_[num];
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::PolygonalTimeSlice::PolygonalTimeSlice(double begin,
                                                   double end,
                                                   double vbegin,
                                                   double vend)
  : TimeSlice(begin,end),
    value_begin_(vbegin),
    value_end_(vend)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::PolygonalTimeSlice::f(double t)
{
  dsassert(contains(t), "wrong time slice called");
  return value_begin_ + (value_end_-value_begin_)/(end()-begin()) * (t-begin());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> DRT::UTILS::PolygonalTimeSlice::FctDer(const double t,
                                                           const unsigned deg)
{
  dsassert(contains(t), "wrong time slice called");

  // resulting vector holding
  std::vector<double> res(deg+1);

  // value of time curve f at time t (0th derivative)
  res[0] = f(t);

  // 1st derivative of time curve f at time t
  if (deg >= 1)
  {
    res[1] = (value_end_-value_begin_)/(end()-begin());
  }

  // 2nd derivative of time curve f at time t
  if (deg >= 2)
  {
    res[2] = 0;
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExplicitTimeSlice::ExplicitTimeSlice(int numex,
                                                 double c1,
                                                 double c2)
  : TimeSlice(0.,1e100),
    numex_(numex),
    c1_(c1),
    c2_(c2)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename ScalarT>
ScalarT DRT::UTILS::ExplicitTimeSlice::Fct(const ScalarT& T)
{
  ScalarT fac = 1;

  switch (numex_)
  {
  case -1: /* f(t)=sin(t:C1*PI:2)_for_t<_C1_else_f(t)=1 */
    if (T <= c1_)
    {
      ScalarT val1 = T/c1_*M_PI/2;
      fac = sin(val1);
    }
    else
      fac = 1.;
    break;
  case -2: /* f(t)=exp(1-1:t)_for_t<C1_else_const. */
    if (T < EPS6)
    {
      fac = 0.;
    }
    else if (T<=c1_ || c1_<EPS6)
    {
      ScalarT val1 = 1. - 1./T;
      fac = exp(val1);
    }
    else
      fac = exp(1. - 1./c1_);
    break;
  case -3: /* f(t)=1-cos(2*PI*C1*t) */
  {
    ScalarT val1 = 2.*M_PI*c1_*T;
    fac  = 1. - cos(val1);
    break;
  }
  case -4: /* f(t)=C2*sin(2PI*C1*t) */
  {
    ScalarT val1 = 2.*c1_*M_PI*T;
    fac  = c2_*sin(val1);
    break;
  }
  case -5: /* f(t)=(sin(PI(t:C1-0.5))+1)*0.5 */
    if (T<=c1_)
    {
      ScalarT val1 = M_PI*(T/c1_-1./2.);
      fac = (sin(val1)+1.)/2.;
    }
    else
      fac = 1.;
    break;
  case -6: /* Beltrami-Flow */
  {
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
    if (id==-1) dserror("Newtonian fluid material could not be found");
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
    ScalarT visc = (ScalarT) actmat->viscosity_;
    ScalarT d = M_PI/2.;
    ScalarT val1 = -c1_*visc*d*d*T;
    fac = exp(val1);
    break;
  }
  case -7: /* Kim-Moin-Flow */
  {
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
    if (id==-1) dserror("Newtonian fluid material could not be found");
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
    ScalarT visc = (ScalarT) actmat->viscosity_;
    ScalarT a = 2.0;
    ScalarT val1 = -c1_*a*a*M_PI*M_PI*visc*T;
    fac = exp(val1);
    break;
  }
  case -8: /* f(t)=(C2/2PI*C1)*cos(2PI*C1*t) +s0*/
  {
    ScalarT val1 = 2.*c1_*M_PI;
    ScalarT s0   = -c2_/val1;
    fac = c2_/val1*cos(val1*T)+s0;
    break;
  }
  case -9: /* f(t)=t:2-C1:(2PI)*cos(PI*t:C1-PI:2) */
    if (T<=c1_)
    {
      ScalarT val1 = M_PI / c1_;
      fac = T*0.5 - 0.5/val1 * cos(val1*T-M_PI*0.5);
    }
    else
      fac = T - c1_ * 0.5;
    break;
  case -11: /* f(t)=-0.5*cos(PI*(T-c1)/c2)+0.5 */
    if (T<c1_)
      fac = 0.0;
    else if ((c1_ <= T) && (T <= (c1_ + c2_)))
      fac = -0.5 * cos(M_PI * ( T - c1_ ) / c2_) + 0.5;
    else
      fac = 1.0;
    break;
  case -12: /* f(t)=(1-c2/c1)*T+c2 (c2: ratio of initial and final value) */
    if (T <= c1_)
      fac = ( ( 1.0 - c2_ ) / c1_ ) * T + c2_;
    else
      fac = 1.0;
    break;
  default:
    dserror("Number of explicit timecurve (NUMEX=%d) unknown", numex_);
  }

  return fac;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::ExplicitTimeSlice::f(double T)
{
  double fac = Fct<double>(T);
  return fac;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> DRT::UTILS::ExplicitTimeSlice::FctDer(const double t,
                                                          const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg+1);

  // function
  res[0] = f(t);

  // derivatives
  if (deg >= 1)
  {
    // number of independent variables
    const int nvar = 1;
    // index of independent variable
    const int ivar = 0;

    // Fad object for evaluation time
    Sacado::Fad::DFad< Sacado::Fad::DFad<double> > tfad(nvar, ivar, t);
    // for 1st order derivative
    tfad.val() = Sacado::Fad::DFad<double>(nvar, ivar, t);

    // 1st & 2nd derivative of time curve function
    Sacado::Fad::DFad<Sacado::Fad::DFad<double> > fdfad
      = Fct< Sacado::Fad::DFad< Sacado::Fad::DFad<double> > >(tfad);

    // return 1st derivative value at time t
    res[1] = fdfad.dx(ivar).val();

    // return 2nd derivative value at time t
    if (deg >= 2)
    {
      res[2] = fdfad.dx(ivar).dx(ivar);
    }
  }

  // deliver function (and its derivatives)
  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::LungTimeSlice::LungTimeSlice(double frequ,
                                         double ppeep,
                                         double phase)
  : TimeSlice(0.,1e100),
    frequ_(frequ),
    ppeep_(ppeep),
    phase_(phase)
{
  if (phase_ == 0.0)
    phase_ = 1/frequ_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename ScalarT>
ScalarT DRT::UTILS::LungTimeSlice::Fct(const ScalarT& t)
{
  if (t <= phase_)
  {
    return 0.5*ppeep_*(1 - cos(1/phase_*M_PI*t));
  }
  else
  {
    return 0.5*(1-ppeep_)*(1 - cos(2*M_PI*frequ_*(t - phase_))) + ppeep_;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::LungTimeSlice::f(double t)
{
  double fac = Fct<double>(t);
  return fac;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> DRT::UTILS::LungTimeSlice::FctDer(const double t,
                                                      const unsigned deg)
{
  //dserror("2nd time differentation is not implemented");

  // resulting vector holding
  std::vector<double> res(deg+1);

  // function
  res[0] = f(t);

  // derivatives
  if (deg >= 1)
  {
    // number of independent variables
    const int nvar = 1;
    // index of independent variable
    const int ivar = 0;

    // Fad object for evaluation time
    Sacado::Fad::DFad< Sacado::Fad::DFad<double> > tfad(nvar, ivar, t);
    // for 1st order derivative
    tfad.val() = Sacado::Fad::DFad<double>(nvar, ivar, t);

    // 1st & 2nd derivative of time curve function
    Sacado::Fad::DFad<Sacado::Fad::DFad<double> > fdfad
      = Fct< Sacado::Fad::DFad< Sacado::Fad::DFad<double> > >(tfad);

    // return 1st derivative value at time t
    res[1] = fdfad.dx(ivar).val();

    // return 2nd derivative value at time t
    if (deg >= 2)
    {
      res[2] = fdfad.dx(ivar).dx(ivar);
    }
  }

  // deliver function (and its derivatives)
  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::BloodTimeSlice::BloodTimeSlice(double period,
                                           double flowrate,
                                           int points,
                                           std::vector<double>& ArrayLength)
  : TimeSlice(0.,1e100),
    period_(period),
    flowrate_(flowrate),
    points_(points),
    ArrayLength_(ArrayLength)

{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename ScalarT>
ScalarT DRT::UTILS::BloodTimeSlice::Fct(const ScalarT& t)
{
  const int DataLength=points_;
  vector<double> SampleNumber;
  vector<double> EvenCoefficient;
  vector<double> OddCoefficient;
  ScalarT fac;
  double C = (double)points_;


  for (int p=0; p<DataLength; p++)
    {
      SampleNumber.push_back(ArrayLength_[p]*flowrate_);
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
        + (ScalarT)EvenCoefficient[h]*cos(2*PI*h*t/period_)
        + (ScalarT)OddCoefficient[h]*sin(2*PI*h*t/period_);
  }


  return fac;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::BloodTimeSlice::f(double t)
{
  double fac = Fct<double>(t);
  return fac;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> DRT::UTILS::BloodTimeSlice::FctDer(const double t,
                                                       const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg+1);

  // function
  res[0] = f(t);

  // derivatives
  if (deg >= 1)
  {
    // number of independent variables
    const int nvar = 1;
    // index of independent variable
    const int ivar = 0;

    // Fad object for evaluation time
    Sacado::Fad::DFad< Sacado::Fad::DFad<double> > tfad(nvar, ivar, t);
    // for 1st order derivative
    tfad.val() = Sacado::Fad::DFad<double>(nvar, ivar, t);

    // 1st & 2nd derivative of time curve function
    Sacado::Fad::DFad<Sacado::Fad::DFad<double> > fdfad
      = Fct< Sacado::Fad::DFad< Sacado::Fad::DFad<double> > >(tfad);

    // return 1st derivative value at time t
    res[1] = fdfad.dx(ivar).val();

    // return 2nd derivative value at time t
    if (deg >= 2)
    {
      res[2] = fdfad.dx(ivar).dx(ivar);
    }
  }

  // deliver function (and its derivatives)
  return res;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprTimeSlice::ExprTimeSlice(double begin, double end, char* buf)
  : TimeSlice(begin,end),
    parsexpr_(DRT::PARSER::Parser<double>(string(buf))),
    parsexprdd_(DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > >(string(buf)))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprTimeSlice::~ExprTimeSlice()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::ExprTimeSlice::f(double t)
{
  dsassert(contains(t), "wrong time slice called");
  //cout << "Curve factor " << parsexpr_.EvaluateCurve(t) << endl;
  return parsexpr_.EvaluateCurve(t);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> DRT::UTILS::ExprTimeSlice::FctDer(const double t,
                                                      const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg+1);

  // function
  res[0] = f(t);
  if (std::isnan(res[0]))
  {
    cout << (*this) << endl;
    cout << "0. derivative at t = "<< t <<": " << res[0] << endl;
    dserror("");
  }

  // derivatives
  if (deg >= 1)
  {
    // number of independent variables
    const int nvar = 1;
    // index of independent variable
    const int ivar = 0;

    // Fad object for evaluation time
    Sacado::Fad::DFad< Sacado::Fad::DFad<double> > tfad(nvar, ivar, t);
    // for 1st order derivative
    tfad.val() = Sacado::Fad::DFad<double>(nvar, ivar, t);

    // 1st & 2nd derivative of time curve function
    Sacado::Fad::DFad<Sacado::Fad::DFad<double> > fdfad
      = parsexprdd_.EvaluateCurve(tfad);

    // return 1st derivative value at time t
    res[1] = fdfad.dx(ivar).val();
    if (std::isnan(res[1]))
    {
      cout << (*this) << endl;
      cout << "1. derivative at t = "<< t <<": " << res[1] << endl;
      dserror("");
    }

    // 2nd derivative requested
    if (deg >= 2)
    {
      // set 2nd derivative at time t
      res[2] = fdfad.dx(ivar).dx(ivar);
      if (std::isnan(res[2]))
      {
        cout << (*this) << endl;
        cout << "2. derivative at t = "<< t <<": " << res[2] << endl;
        dserror("");
      }
    }
  }

  // return function (and its derivatives)
  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::TimeCurve::f(double t)
{
  if (slices_.size()==0)
    dserror("No time slices defined. Fix input.");

  // if the requested time is before the start we use the first
  // available time we have
  RefCountPtr<TimeSlice> slice = slices_[0];
  if (t < slice->begin())
    return slice->f(slice->begin());

  // look for the right slice and ask it
  for (unsigned i=0; i<slices_.size(); ++i)
  {
    slice = slices_[i];
    if (slice->contains(t))
      return slice->f(t);
    if (slice->begin() > t)
      dserror("a gap between time slices occured");
  }

  // if we exceed our slices use the last available time
  slice = slices_.back();
  return slice->f(slice->end());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> DRT::UTILS::TimeCurve::FctDer(const double t,
                                                  const unsigned deg)
{
  // verify
  if (deg > 2)
    dserror("Highest degree of differentation must be in [0,2]");

  if (slices_.size()==0)
    dserror("No time slices defined. Fix input.");

  // if the requested time is before the start we use the first
  // available time we have
  RefCountPtr<TimeSlice> slice = slices_[0];
  if (t < slice->begin())
  {
    return slice->FctDer(slice->begin(), deg);
  }

  // look for the right slice and ask it
  for (unsigned i=0; i<slices_.size(); ++i)
  {
    slice = slices_[i];
    if (slice->contains(t))
    {
      return slice->FctDer(t, deg);
    }
    if (slice->begin() > t)
      dserror("a gap between time slices occured");
  }

  // if we exceed our slices use the last available time
  slice = slices_.back();
  return slice->FctDer(slice->end(), deg);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::TimeCurve::end()
{
  return slices_[slices_.size()-1]->end();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::TimeCurve::AddSlice(Teuchos::RefCountPtr<TimeSlice> slice)
{
  // Do we need more error checking here?
  slices_.push_back(slice);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::PolygonalTimeSlice::Print(std::ostream& out) const
{
  out << "    PolygonalTimeSlice(begin=" << begin()
      << ", end=" << end()
      << ", value_begin=" << value_begin_
      << ", value_end=" << value_end_
      << ")\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::ExplicitTimeSlice::Print(std::ostream& out) const
{
  out << "    ExplicitTimeSlice(numex=" << numex_
      << ", c1=" << c1_
      << ", c2=" << c2_
      << ")\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::LungTimeSlice::Print(std::ostream& out) const
{
  out << "    LungTimeSlice(frequ=" << frequ_
      << ", ppeep=" << ppeep_
      << ", phase=" << phase_
      << ")\n";
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void DRT::UTILS::BloodTimeSlice::Print(std::ostream& out) const
{
  out << "   BloodTimeSlice(period=" << period_
      << ",  flowrate=" << flowrate_
      << ")\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::ExprTimeSlice::Print(std::ostream& out) const
{
  out << "    ExprTimeSlice(begin=" << begin()
      << ", end=" << end()
      << ")\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& DRT::UTILS::operator<<(std::ostream& out, const DRT::UTILS::TimeSlice& slice)
{
  slice.Print(out);
  return out;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& DRT::UTILS::operator<<(std::ostream& out, const DRT::UTILS::TimeCurve& curve)
{
  out << "  Time Curve:\n";
  for (unsigned i=0; i<curve.slices_.size(); ++i)
  {
    out << *curve.slices_[i];
  }
  return out;
}

#endif
