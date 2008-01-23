/*----------------------------------------------------------------------*/
/*!
\file drt_timecurve.cpp

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

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>

#include "drt_timecurve.H"
#include "drt_dserror.H"


#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


using namespace std;
using namespace Teuchos;


/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*/
// the static instance
/*----------------------------------------------------------------------*/
DRT::UTILS::TimeCurveManager DRT::UTILS::TimeCurveManager::instance_;


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
    if (frfind(const_cast<char*>(curve.str().c_str()))==1)
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
          char* colpointer = strstr(allfiles.actplace,"BYABSTIME");
          colpointer += 13;
          double begin = strtod(colpointer,&colpointer);
          double end = strtod(colpointer,&colpointer);
          colpointer = strstr(allfiles.actplace,"FACTOR");
          colpointer += 6;
          double vbegin = strtod(colpointer,&colpointer);
          double vend = strtod(colpointer,&colpointer);
          curve.AddSlice(rcp(new PolygonalTimeSlice(begin,end,vbegin,vend)));
          continue;
        }
        frchk("Explicit",&ierr);
        if (ierr==1)
        {
          char buffer[50];
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
          else
            dserror("cannot read function of CURVE%d",i);

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
          char buffer[250];
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
        frchk("PhysiologicalWavefrom",&ierr);
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
DRT::UTILS::ExplicitTimeSlice::ExplicitTimeSlice(int numex, double c1, double c2)
  : TimeSlice(0.,1e100),
    numex_(numex),
    c1_(c1),
    c2_(c2)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::ExplicitTimeSlice::f(double T)
{
  double fac = 1.0;

  switch (numex_)
  {
  case -1: /* f(t)=sin(t:C1*PI:2)_for_t<_C1_else_f(t)=1 */
    if (T <= c1_)
    {
      double val1 = T/c1_*M_PI/2;
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
      double val1 = 1. - 1./T;
      fac = exp(val1);
    }
    else
      fac = exp(1. - 1./c1_);
    break;
  case -3: /* f(t)=1-cos(2*PI*C1*t) */
  {
    double val1 = 2.*M_PI*c1_*T;
    fac  = 1. - cos(val1);
    break;
  }
  case -4: /* f(t)=C2*sin(2PI*C1*t) */
  {
    double val1 = 2.*c1_*M_PI*T;
    fac  = c2_*sin(val1);
    break;
  }
  case -5: /* f(t)=(sin(PI(t:C1-0.5))+1)*0.5 */
    if (T<=c1_)
    {
      double val1 = M_PI*(T/c1_-1./2.);
      fac = (sin(val1)+1.)/2.;
    }
    else
      fac = 1.;
    break;
  case -6: /* Beltrami-Flow */
  {
    double visc = mat[0].m.fluid->viscosity;
    double d = M_PI/2.;
    double val1 = -c1_*visc*d*d*T;
    fac = exp(val1);
    break;
  }
  case -7: /* Kim-Moin-Flow */
  {
    double visc = mat[0].m.fluid->viscosity;
    double a = 2.0;
    double val1 = -c1_*a*a*M_PI*M_PI*visc*T;
    fac = exp(val1);
    break;
  }
  case -8: /* f(t)=(C2/2PI*C1)*cos(2PI*C1*t) +s0*/
  {
    double val1 = 2.*c1_*M_PI;
    double s0   = -c2_/val1;
    fac = c2_/val1*cos(val1*T)+s0;
    break;
  }
  case -9: /* f(t)=t:2-C1:(2PI)*cos(PI*t:C1-PI:2) */
    if (T<=c1_)
    {
      double val1 = M_PI / c1_;
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
  default:
    dserror("Number of explicit timecurve (NUMEX=%d) unknown", numex_);
  }

  return fac;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::LungTimeSlice::LungTimeSlice(double frequ, double ppeep, double phase)
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
double DRT::UTILS::LungTimeSlice::f(double t)
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
DRT::UTILS::BloodTimeSlice::BloodTimeSlice(double period, double flowrate, int points,  std::vector<double>& ArrayLength )
  : TimeSlice(0.,1e100),
    period_(period),
    flowrate_(flowrate),
    points_(points),
    ArrayLength_(ArrayLength)

{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::BloodTimeSlice::f(double t)
{


  const int DataLength=points_;
 //  double EvenCoefficient[DataLength/2]={0};
//   double OddCoefficient[DataLength/2]={0};
//   double SampleNumber[DataLength]={0};
   double EvenCoefficient[31]={0};
  double OddCoefficient[31]={0};
  double SampleNumber[60]={0};
  double fac;
  double C = (double)points_;

  // printf("%d\n",DataLength);

  for (int p=0;p<DataLength;p++){
    SampleNumber[p]=ArrayLength_[p]*flowrate_;
      }

  for (int p=0; p<=DataLength/2; p++){
   EvenCoefficient[p] = 0;
   OddCoefficient[p] = 0;

   for (int num=0; num<=DataLength-1; num++){
     EvenCoefficient[p] = EvenCoefficient[p]+2/C*SampleNumber[num]*cos(2*PI*p*(num+1)/C);
     OddCoefficient[p] = OddCoefficient[p]+2/C*SampleNumber[num]*sin(2*PI*p*(num+1)/C);
   }
   //  printf("%3d : % f  % f\n", p, EvenCoefficient[p], OddCoefficient[p]);
 }

EvenCoefficient[DataLength/2] = EvenCoefficient[DataLength/2]/2;
OddCoefficient[DataLength/2] = 0;
fac = EvenCoefficient[0]/2;

  for (int h=1; h<=DataLength/2; h++){
    fac = fac+EvenCoefficient[h]*cos(2*PI*h*t/period_)+OddCoefficient[h]*sin(2*PI*h*t/period_);
  }

  return fac;
}




/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprTimeSlice::ExprTimeSlice(double begin, double end, char* buf)
  : TimeSlice(begin,end),
    expr_(pss_parse(buf))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprTimeSlice::~ExprTimeSlice()
{
  pss_parse_cleanup(expr_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::ExprTimeSlice::f(double t)
{
  dsassert(contains(t), "wrong time slice called");
  return pss_evaluate_curve(expr_,t);
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& DRT::UTILS::operator<<(std::ostream& out, const DRT::UTILS::TimeCurveManager& manager)
{
  out << "Time Curve Manager:\n";
  for (unsigned i=0; i<manager.curves_.size(); ++i)
  {
    out << manager.curves_[i];
  }
  return out;
}

#endif
