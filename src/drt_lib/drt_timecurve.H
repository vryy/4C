/*----------------------------------------------------------------------*/
/*!
\file drt_timecurve.H

\brief Managing and evaluating of time curves

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef DRT_TIMECURVE_H
#define DRT_TIMECURVE_H

#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>

#include "drt_dserror.H"
#include "drt_parser.H"

extern "C"
{
#include "../pss_full/pss_parser.h"
}

namespace DRT {
namespace UTILS {


  /// base of all slices of a time curve with start and end time
  class TimeSlice
  {
  public:

    /// construct with begin and end time
    TimeSlice(double begin, double end) : begin_(begin), end_(end) {}

    /// virtual destructor because this class will be derived
    virtual ~TimeSlice() {}

    /// interval check [begin:end[
    bool contains(double t) const { return t>=begin_-1e-13 and t<end_+1e-13; }

    /// evaluate time curve at given time
    virtual double f(double time) = 0;

    /// begin time of slice
    double begin() const { return begin_; }

    /// end time of slice
    double end() const { return end_; }

    /// debug output of this slice
    virtual void Print(std::ostream& out) const = 0;

  private:
    double begin_;
    double end_;

    /// output operator
    friend std::ostream& operator<<(std::ostream& out, const TimeSlice& slice);
  };


  /// time slice defined by polygonal (with jumps!)
  class PolygonalTimeSlice : public TimeSlice
  {
  public:
    PolygonalTimeSlice(double begin, double end, double vbegin, double vend);

    /// evaluate time curve at given time
    double f(double t);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;
  private:

    /// value at begin of time slice
    double value_begin_;

    /// value at end of time slice
    double value_end_;
  };



  /// explicit strings with fixed implementation
  /*!
    Depreciated. Do not extend this class. Create a new one if you
    need. Or even better use ExprTimeSlice.
  */
  class ExplicitTimeSlice : public TimeSlice
  {
  public:
    ExplicitTimeSlice(int numex, double c1, double c2);

    /// evaluate time curve at given time
    double f(double t);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;
  private:
    int numex_;
    double c1_;
    double c2_;
  };



  /// special respiration time slice
  class LungTimeSlice : public TimeSlice
  {
  public:
    LungTimeSlice(double frequ, double ppeep, double phase);

    /// evaluate time curve at given time
    double f(double t);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;
  private:
    double frequ_;
    double ppeep_;
    double phase_;
  };

   /// special respiration time slice
  class BloodTimeSlice : public TimeSlice
  {
  public:
    BloodTimeSlice(double period, double flowrate, int points,  std::vector<double>& ArrayLength );

    /// evaluate time curve at given time
    double f(double t);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;
  private:
    double period_;
    double flowrate_;
    int points_;
    std::vector<double> ArrayLength_;
  };





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

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;

  private:

    /// expression syntax tree
    _ST_NODE* expr_;

    /// parsed expression as syntax tree
    DRT::Parser<double> parexpr_;
  };



  /// time curve consisting of a bunch of time slices
  /*!
    A time curve consisting of one or more (adjacent) time
    slices. Each time slice can be of a different type, so you can
    have a total mess of a curve. But honestly, the slice class you
    will normally want to use is ExprTimeSlice. The other ones are
    provided to stay compatible with old code only.

    All time curves are managed by TimeCurveManager.
   */
  class TimeCurve
  {
    friend class TimeCurveManager;
  public:

    /// evaluate time curve at given time
    /*!
      If time t is less than the starting time of the first slice,
      this starting time is used instead. If time t is larger than the
      ending time of the last slice, this ending time is used.
     */
    double f(double t);

    // return the end time of the last time slice of this curve
    double end() {return slices_[slices_.size()-1]->end();}

  private:

    /// add a new time slice at the and of this curve
    void AddSlice(Teuchos::RefCountPtr<TimeSlice> slice);

    /// all our time slices ordered by time
    std::vector<Teuchos::RefCountPtr<TimeSlice> > slices_;

    /// output operator
    friend std::ostream& operator<<(std::ostream& out, const TimeCurve& curve);
  };



  /// Singleton class that manages all defined time curves
  /*!
    There is one instance of this class and it contains all the time
    curves that have been defined. You can get this instance via the
    static Instance() method and ask it for the curve you are
    interessted it. So what we have here is a global variable in
    disguise. After all, that's why it is called a Singleton class.

    Time curves consist of one or more slices where each slice has a
    specified starting time and ending time. This way you can specify
    discontinous time curves easily. Look at TimeCurve for further
    information.
   */
  class TimeCurveManager
  {
  public:

    /// return the single instance we have of this class
    static TimeCurveManager& Instance() { return instance_; }

    /// read the baci input file and set up all curves
    void ReadInput();

    /// return a specified time curve
    TimeCurve& Curve(int num) { return curves_[num]; }

  private:

    /// private default constructor to disallow creation of instances
    TimeCurveManager() {}

    /// no copying
    TimeCurveManager(const TimeCurveManager&);

    /// no copying
    TimeCurveManager& operator=(const TimeCurveManager&);

    /// all curves defined
    std::vector<TimeCurve> curves_;

    /// the single instance
    static TimeCurveManager instance_;

    /// output operator
    friend std::ostream& operator<<(std::ostream& out, const TimeCurveManager& manager);
  };

}  // end namespace UTILS
}  // end namespace DRT

#endif  // DRT_TIMECURVE_H
#endif  // CCADISCRET
