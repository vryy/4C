/*----------------------------------------------------------------------*/
/*!
\file drt_timecurve.cpp

\brief Implementation Managing and evaluating of time curves
<pre>
\brief Implementation
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <Sacado.hpp>
#include <cmath>

#include "drt_parser.H"
#include "drt_timecurve.H"
#include "drt_dserror.H"
#include "drt_linedefinition.H"

#include "drt_globalproblem.H"
#include "standardtypes_cpp.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matpar_bundle.H"


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

    /// \brief Implementation Evaluate time curve and its derivatives
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

  std::ostream& operator<<(std::ostream& out, const TimeSlice& slice);
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
                   std::string type_of_interp,
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
    std::string interp_type_;
  };


  /*--------------------------------------------------------------------*/
  /// special time slice for periodic repetition of other time slices
  class PeriodicTimeSlice : public TimeSlice
  {
  public:
    PeriodicTimeSlice(
        const double begin,
        const double end,
        const double period,
        DRT::UTILS::TimeCurve& curve);

    /// evaluate time curve at given time
    double f(double t);

    /// evaluate reference time in the first periodic cycle
    double GetReferenceTime(double t);

    /// evaluate time curve and its derivatives
    std::vector<double> FctDer(const double t, const unsigned deg);

    /// debug output of this slice
    virtual void Print(std::ostream& out) const;

  private:
    // member variables begin_ and end_ are inherited from TimeSlice !!
    const double  period_; // length of period
    TimeCurve&    curve_;  // reference to time curve that owns this time slice
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
    There might be applications that require this flexibility.
   */
  class ExprTimeSlice : public TimeSlice
  {
  public:
    /// construct syntax tree from std::string buffer
    ExprTimeSlice(double begin, double end, std::string buf);

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
//! Print function
/*----------------------------------------------------------------------*/
void PrintTimeCurveDatHeader()
{
  DRT::UTILS::TimeCurveManager timecurvemanager;
  Teuchos::RCP<DRT::INPUT::Lines> lines = timecurvemanager.ValidTimeCurveLines();

  lines->Print(std::cout);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::INPUT::Lines> DRT::UTILS::TimeCurveManager::ValidTimeCurveLines()
{
  DRT::INPUT::LineDefinition polygonal;
  polygonal
    .AddNamedInt("CURVE")

    // This is a really stupid definition. Once upon a time there were
    // different options, but since those are no longer supported, we do not
    // care.
    .AddTag("on")
    .AddTag("Polygonal")
    .AddNamedDouble("T")
    .AddTag("BYSTEP")
    .AddNamedDoubleVector("No",2)
    .AddTag("BYABSTIME")
    .AddNamedDoubleVector("Yes",2)
    .AddNamedDoubleVector("FACTOR",2)
    ;

  DRT::INPUT::LineDefinition smallpolygonal;
  smallpolygonal
    .AddNamedInt("CURVE")

    // This is a really stupid definition. Once upon a time there were
    // different options, but since those are no longer supported, we do not
    // care.
    .AddTag("on")
    .AddTag("Polygonal")
    .AddNamedDouble("T")
    .AddTag("BYABSTIME")
    .AddNamedDoubleVector("Yes",2)
    .AddNamedDoubleVector("FACTOR",2)
    ;

  DRT::INPUT::LineDefinition expl;
  expl
    .AddNamedInt("CURVE")
    .AddTag("on")
    .AddTag("Explicit")
    .AddNamedString("FUNC")
    .AddNamedDouble("c1")
    .AddNamedDouble("c2")
    ;

  DRT::INPUT::LineDefinition func;
  func
    .AddNamedInt("CURVE")
    .AddTag("on")
    .AddTag("EXPR")
    .AddNamedString("FUNC")
    .AddNamedDouble("t1")
    .AddNamedDouble("t2")
    ;

  DRT::INPUT::LineDefinition lungsinus;
  lungsinus
    .AddNamedInt("CURVE")
    .AddTag("on")
    .AddTag("LungSinus")
    .AddNamedDouble("Frequ")
    .AddNamedDouble("pPEEP")
    .AddNamedDouble("Phase")
    ;

  DRT::INPUT::LineDefinition physiologicalwaveform;
  physiologicalwaveform
    .AddNamedInt("CURVE")
    .AddTag("on")
    .AddTag("PhysiologicalWaveform")
    .AddNamedString("InterpolType")
    .AddNamedDouble("Period")
    .AddNamedDouble("Flowrate")
    .AddNamedInt("Samplingpoints")
    .AddNamedDoubleVector("Arrayread","Samplingpoints")
    ;

  DRT::INPUT::LineDefinition periodicrepetition;
  periodicrepetition
    .AddNamedInt("CURVE")
    .AddTag("on")
    .AddTag("PeriodicRepetition")
    .AddNamedDouble("Period")
    .AddNamedDouble("t1")
    .AddNamedDouble("t2")
    ;

  DRT::INPUT::LineDefinition curvefromfile;
  curvefromfile
    .AddNamedInt("CURVE")
    .AddTag("on")
    .AddTag("File")
    .AddNamedString("PathToFile")
    ;

  Teuchos::RCP<DRT::INPUT::Lines> lines = Teuchos::rcp(new DRT::INPUT::Lines("CURVE"));
  lines->Add(polygonal);
  lines->Add(smallpolygonal);
  lines->Add(expl);
  lines->Add(func);
  lines->Add(lungsinus);
  lines->Add(physiologicalwaveform);
  lines->Add(periodicrepetition);
  lines->Add(curvefromfile);
  return lines;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::TimeCurveManager::ReadInput(DRT::INPUT::DatFileReader& reader)
{
  curves_.clear();

  Teuchos::RCP<DRT::INPUT::Lines> lines = ValidTimeCurveLines();

  // test for as many curves as there are
  for (int i = 1;; ++i)
  {
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition> > curves = lines->Read(reader, i);
    if (curves.size()==0)
      break;

    // so we have a new time curve
    curves_.push_back(Teuchos::rcp(new TimeCurve()));
    TimeCurve& curve = *(curves_.back());

    for (unsigned j=0; j<curves.size(); ++j)
    {
      int id;
      curves[j]->ExtractInt("CURVE",id);
      if (id!=i) dserror("expected CURVE%d but got CURVE%d", i, id);

      if (curves[j]->HaveNamed("Polygonal"))
      {
        std::vector<double> byabstime;
        curves[j]->ExtractDoubleVector("Yes",byabstime);
        std::vector<double> factor;
        curves[j]->ExtractDoubleVector("FACTOR",factor);
        curve.AddSlice(Teuchos::rcp(new PolygonalTimeSlice(byabstime[0],byabstime[1],
                                                  factor[0],factor[1])));
      }
      else if (curves[j]->HaveNamed("Explicit"))
      {
        std::string buffer;
        curves[j]->ExtractString("FUNC",buffer);

        int numex = 0;
        if (std::string(buffer)=="f(t)=sin(t:C1*PI:2)_for_t<_C1_else_f(t)=1")
          numex=-1;
        else if (std::string(buffer)=="f(t)=exp(1-1:t)_for_t<C1_else_const.")
          numex=-2;
        else if (std::string(buffer)=="f(t)=1-cos(2*PI*C1*t)")
          numex=-3;
        else if (std::string(buffer)=="f(t)=C2*sin(2PI*C1*t)")
          numex=-4;
        else if (std::string(buffer)=="f(t)=(sin(PI(t:C1-0.5))+1)*0.5")
          numex=-5;
        else if (std::string(buffer)=="BELTRAMI")
          numex=-6;
        else if (std::string(buffer)=="KIM-MOIN")
          numex=-7;
        else if (std::string(buffer)=="f(t)=(C2:2PI*C1)*cos(2PI*C1*t)")
          numex=-8;
        else if (std::string(buffer)=="f(t)=t:2-C1:(2PI)*cos(PI*t:C1-PI:2)")
          numex=-9; /* time integral of numex -5 */
        else if (std::string(buffer)=="f(t)=-0.5*cos(PI*(T-C1)/C2)+0.5")
          numex=-11; /* ramp function with horizontal slopes */
        else if (std::string(buffer)=="f(t)=(1-C2/C1)*T+C2")
          numex=-12; /* linearly increasing function with non-zero initial value */
        else if (std::string(buffer)=="f(t)=0.5+0.5*cos(PI*(T-C1)/(C2-C1))")
          numex=-13;
        else
          dserror("Cannot read function of CURVE%d: %s",i,std::string(buffer).c_str());

        double c1;
        double c2;
        curves[j]->ExtractDouble("c1",c1);
        curves[j]->ExtractDouble("c2",c2);

        curve.AddSlice(Teuchos::rcp(new ExplicitTimeSlice(numex,c1,c2)));
      }
      else if (curves[j]->HaveNamed("EXPR"))
      {
        std::string buffer;
        curves[j]->ExtractString("FUNC",buffer);

        double begin;
        double end;

        curves[j]->ExtractDouble("t1",begin);
        curves[j]->ExtractDouble("t2",end);

        curve.AddSlice(Teuchos::rcp(new ExprTimeSlice(begin, end, buffer)));
      }
      else if (curves[j]->HaveNamed("LungSinus"))
      {
        double frequ;
        double ppeep;
        double phase;

        curves[j]->ExtractDouble("Frequ",frequ);
        curves[j]->ExtractDouble("pPEEP",ppeep);
        curves[j]->ExtractDouble("Phase",phase);

        curve.AddSlice(Teuchos::rcp(new LungTimeSlice(frequ, ppeep, phase)));
      }
      else if (curves[j]->HaveNamed("PhysiologicalWaveform"))
      {
        double period;
        double flowrate;
        int points;
        std::vector<double> Arrayread;
        std::string type_of_interp;

        curves[j]->ExtractString("InterpolType",type_of_interp);
        curves[j]->ExtractDouble("Period",period);
        curves[j]->ExtractDouble("Flowrate",flowrate);
        curves[j]->ExtractInt("Samplingpoints",points);
        curves[j]->ExtractDoubleVector("Arrayread",Arrayread);

        curve.AddSlice(Teuchos::rcp(new BloodTimeSlice(period, flowrate,type_of_interp, points, Arrayread)));
      }
      else if (curves[j]->HaveNamed("PeriodicRepetition"))
      {
        double period;
        double begin;
        double end;

        curves[j]->ExtractDouble("Period",period);
        curves[j]->ExtractDouble("t1",begin);
        curves[j]->ExtractDouble("t2",end);

        curve.AddSlice(Teuchos::rcp(new PeriodicTimeSlice(begin,end,period,curve)));
      }
      else if (curves[j]->HaveNamed("File"))
      {
        std::string path;
        curves[j]->ExtractString("PathToFile",path);

        // check existence of filename
        if (path.size() <= 0)
          dserror("path to time curve file invalid");

        std::string contents;
        std::vector<char> data;
        // read file from disk on proc 0 and broadcast it to all other procs
        if(reader.Comm()->MyPID() == 0)
        {
          // make path relative to input file path
          if (path[0]!='/')
          {
            std::string filename = reader.MyInputfileName();
            std::string::size_type pos = filename.rfind('/');
            if (pos!=std::string::npos)
            {
              std::string tmp = filename.substr(0,pos+1);
              path.insert(path.begin(), tmp.begin(), tmp.end());
            }
          }

          {
            std::ifstream infile(path.c_str());
            if(not infile.is_open())
              dserror("time curve file could not be opened");

            infile.seekg(0, std::ios::end);
            contents.resize(infile.tellg());
            infile.seekg(0, std::ios::beg);
            infile.read(&contents[0], contents.size());
            infile.close();
          }

          // broadcast
          if (reader.Comm()->NumProc() > 1)
          {
          DRT::PackBuffer buffer;
          DRT::ParObject::AddtoPack(buffer, contents);
          buffer.StartPacking();
          DRT::ParObject::AddtoPack(buffer, contents);
          std::swap(data, buffer());
          }
        }

        // broadcast
        if (reader.Comm()->NumProc() > 1)
        {
          ssize_t data_size = data.size();
          reader.Comm()->Broadcast(&data_size,1,0);
          if (reader.Comm()->MyPID() != 0)
            data.resize(data_size,0);
          reader.Comm()->Broadcast(&(data[0]), data.size(), 0);

          if (reader.Comm()->MyPID() != 0)
          {
            size_t pos = 0;
            DRT::ParObject::ExtractfromPack(pos, data, contents);
          }
        }
        data.clear();

        // start processing the input file
        double t0;
        double t1;
        double val0;
        double val1;

        std::string line;
        std::istringstream iss(contents);
        // read until first valid line
        while(getline(iss, line))
        {
          std::stringstream ss0(line);
          ss0 >> t0 >> val0;
          if(not ss0.fail())
            break;
        }

        if(iss.eof())
          dserror("could not find at least two valid lines in file of CURVE%d", i);

        // read out file
        while(getline(iss, line))
        {
          std::stringstream ss1(line);
          ss1 >> t1 >> val1;

          // skip invalid lines (hopefully comments)
          if(ss1.fail())
            continue;

          // add polygonial time slice
          curve.AddSlice(Teuchos::rcp(new PolygonalTimeSlice(t0,t1,val0,val1)));

          t0 = t1;
          val0 = val1;
        }

      }
      else
        dserror("unknown type of time curve in CURVE%d", i);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::TimeCurve& DRT::UTILS::TimeCurveManager::Curve(int num)
{
  // ensure that desired curve is available (prevents segmentation fault)
  if ((curves_.size()< (unsigned int)(num+1)) || num<0)
    dserror("time curve %d not available \n"
        "possible problems: \n"
        "a) Dirichlet or Neumann condition: \n "
        "   defined NUMDOF's does not match number of dof's defined by "
        "the problem type / material / space dimensions \n"
        "b) The numbering of the curves have to start with 1 ",num+1);

  return *(curves_[num]);
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
    if (c1_<EPS13) dserror ("Illegal constant C1 in time curve");
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
    if (c1_<EPS13) dserror ("Illegal constant C1 in time curve");
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
    // get kinematic viscosity
    ScalarT visc = (ScalarT) actmat->viscosity_ / actmat->density_;
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
    // get kinematic viscosity
    ScalarT visc = (ScalarT) actmat->viscosity_ / actmat->density_;
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
  case -13: /* f(t)=0.5+0.5*cos(PI*(T-C1)/(C2-C1)) */
    if (T<= c1_)
      fac = 1.0;
    else if (T>= c2_)
      fac = 0.0;
    else
      fac = 0.5 + 0.5*cos(M_PI*(T-c1_)/(c2_-c1_));
    break;
  default:
    dserror("Number of explicit timecurve (NUMEX=%d) unknown", numex_);
    break;
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
                                           std::string interp_type,
                                           int points,
                                           std::vector<double>& ArrayLength)
  : TimeSlice(0.,1e100),
    period_(period),
    flowrate_(flowrate),
    points_(points),
    ArrayLength_(ArrayLength),
    interp_type_(interp_type)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <typename ScalarT>
ScalarT DRT::UTILS::BloodTimeSlice::Fct(const ScalarT& t)
{
  const int DataLength=points_;
  std::vector<double> SampleNumber;
  std::vector<double> EvenCoefficient;
  std::vector<double> OddCoefficient;
  ScalarT fac = 0.0;
  double C = (double)points_;

  if (interp_type_ == "Fourier")
  {
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
  }
  else if (interp_type_ == "Linear")
  {
    //double phi = PI/11.0;
    //return -500.0 -350.0/2.0*(1.0 - cos(2.0*PI*t/period_ + phi ));
    // get the time within a period
    double loc_time = Sacado::Fad::DFad<Sacado::Fad::DFad<double> > (t) .val().val();

    if (loc_time < 0.0)
    {
      loc_time += period_;
    }

    loc_time = fmod(loc_time,period_);

    // get the time step size between to array values
    double dt = period_/((double)DataLength - 1.0);

    int index_1   = (int)floor(loc_time/dt);
    double time_1 = double(index_1)*dt;
    double time_2  = time_1 + dt;

    int index_2;

    // get values the interpolation values
    double val1 = ArrayLength_[index_1];

    if (index_1 != DataLength - 1)
    {
      index_2 = index_1 + 1;
    }
    else
    {
      index_2 = 1;
    }
    double val2 = ArrayLength_[index_2];

    //    cout<<"time: "<<t<<" loc_time "<<loc_time<<" time1: "<<time_1<<" time2: "<< time_2<<endl;
    //    cout<<"index_1: "<<index_1<<" index_2: "<<index_2<<" val1: "<<val1<<" val_2: "<<val2<<endl;
    // calculate the interpolated value
    fac  =  (time_2 - loc_time)*val1 + (loc_time - time_1)*val2;
    fac /=  (time_2 - time_1);
    fac *= flowrate_;

  }
  else
  {
    dserror("[%s]: \"PhysiologicalWaveform\" can have only \"Fourier\" or \"Linear\" interpolations",interp_type_.c_str());
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
DRT::UTILS::PeriodicTimeSlice::PeriodicTimeSlice(
    const double begin,
    const double end,
    const double period,
    DRT::UTILS::TimeCurve& curve)
: TimeSlice(begin,end),
    period_(period),
    curve_(curve)
{
  // safety check
  if (period_<EPS13) dserror("Periodic length of size zero or negative");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::PeriodicTimeSlice::GetReferenceTime(double t)
{
  double reftime(t);
  // we go back to the time slices that we want to repeat
  while (reftime > begin())
  {
    reftime -= period_;
  }
  return reftime; // return shifted time value
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::PeriodicTimeSlice::f(double t)
{
  return curve_.f(GetReferenceTime(t));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> DRT::UTILS::PeriodicTimeSlice::FctDer(const double t,
                                                       const unsigned deg)
{
  return curve_.FctDer(GetReferenceTime(t),deg);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprTimeSlice::ExprTimeSlice(double begin, double end, std::string buf)
  : TimeSlice(begin,end),
    parsexpr_(DRT::PARSER::Parser<double>(buf)),
    parsexprdd_(DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > >(buf))
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
    std::cout << (*this) << std::endl;
    std::cout << "0. derivative at t = "<< t <<": " << res[0] << std::endl;
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
      std::cout << (*this) << std::endl;
      std::cout << "1. derivative at t = "<< t <<": " << res[1] << std::endl;
      dserror("");
    }

    // 2nd derivative requested
    if (deg >= 2)
    {
      // set 2nd derivative at time t
      res[2] = fdfad.dx(ivar).dx(ivar);
      if (std::isnan(res[2]))
      {
        std::cout << (*this) << std::endl;
        std::cout << "2. derivative at t = "<< t <<": " << res[2] << std::endl;
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
  Teuchos::RCP<TimeSlice> slice = slices_[0];
  if (t < slice->begin())
    return slice->f(slice->begin());

  // look for the right slice and ask it
  for (unsigned i=0; i<slices_.size(); ++i)
  {
    slice = slices_[i];
    if (slice->contains(t))
      return slice->f(t);
    if (slice->begin() > t)
      dserror("a gap between time slices occurred");
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
  Teuchos::RCP<TimeSlice> slice = slices_[0];
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

  // if we exceed our slices we use the last available time slice
  // for evaluation of f(t). This way the curve becomes a constant function
  // after the last slice. Consequently, first and second derivatives
  // have to be zero in that case to be consistent with f(t) evaluation!
  slice = slices_.back();
  std::vector<double> res = slice->FctDer(slice->end(), deg);
  // set derivatives to zero while f(t) stored in res[0] is kept
  for (size_t i=1; i < res.size();i++)
  {
    res[i] = 0.0;
  }
  return res;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::TimeCurve::end()
{
  return slices_[slices_.size()-1]->end();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::TimeCurve::AddSlice(Teuchos::RCP<TimeSlice> slice)
{
  // Do we need more error checking here?
  slices_.push_back(slice);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::TimeCurve::Print(std::ostream& out) const
{
  out << "  Time Curve:\n";
  for (unsigned i=0; i<slices_.size(); ++i)
  {
    out << *slices_[i];
  }
  return;
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
void DRT::UTILS::PeriodicTimeSlice::Print(std::ostream& out) const
{
  out << "    PeriodicTimeSlice(begin=" << begin()
      << ", end=" << end()
      << ", period=" << period_
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
std::ostream& operator<<(std::ostream& out, const DRT::UTILS::TimeCurve& curve)
{
  curve.Print(out);
  return out;
}

