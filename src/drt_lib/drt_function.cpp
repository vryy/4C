/*----------------------------------------------------------------------*/
/*!
\file drt_function.cpp

\brief Managing and evaluating of spatial functions

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

#include "drt_discret.H"
#include "drt_function.H"
#include "drt_globalproblem.H"
#include "drt_timecurve.H"
#include "drt_linedefinition.H"
#include "../drt_mat/newtonianfluid.H"


namespace DRT {
namespace UTILS {

  /// spatial function based on parsed expression
  class ExprFunction : public Function
  {
  public:
    /*!

    \brief construct spatial function from expression with given origin

    \note  Upon construction, the object defines a spatial function
           returning the same function value for every spatial dimension.
           If a vector-valued spatial function is required, further
           expressions can be added via the AddExpr function. In this
           case, the spatial function defined upon construction will be
           the first component of the vector-valued spatial function.

    \param buf (i) (c-string) expression to be parsed during evaluation
    \param x   (i) x-coordinate of the origin of the coordinate system
    \param y   (i) y-coordinate of the origin of the coordinate system
    \param z   (i) z-coordinate of the origin of the coordinate system

    */
    ExprFunction(char* buf, double x, double y, double z);

    /*!

    \brief Default constructor creating empty object. Expressions are
           added with add function

    */
    ExprFunction();

    /*!

    \brief clean up parse tree

    */
    ~ExprFunction();


    /*!

    \brief evaluate function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functionsb, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief add expression to an existing ExprFunction in order to extend
           it to a vector-valued spatial function.

           Every call to AddExpr adds one more component to the
           vector-valued function.

    \param buf (i) (c-string) expression to be parsed during evaluation of this component
    \param x   (i) x-coordinate of the origin of the coordinate system of this component
    \param y   (i) y-coordinate of the origin of the coordinate system of this component
    \param z   (i) z-coordinate of the origin of the coordinate system of this component

    */
    void AddExpr(std::string buf, double x, double y, double z);
    /*!

    \brief Return the number of components of this spatial function
    (1 for scalar functions, dim for vector-valued functions)

    \return number of components

    */
    virtual int NumberComponents()
      {
        return(expr_.size());
      };

  private:

    /*
      for scalar spatial functions returning the same value for all
      dimensions:

     -----------------------------------------+-------+-------+-------
      spatial function for dimension          |   0   |   1   |   2
     -----------------------------------------+-------+-------+-------
      origin of spatial function, x-component |         x_[0]
     -----------------------------------------+-----------------------
      origin of spatial function, y-component |         y_[0]
     -----------------------------------------+-------+-------+-------
      origin of spatial function, x-component |         z_[0]
     -----------------------------------------+-------+-------+-------


      for vector-valued spatial functions, returning seperate values
      for all dimensions:

     -----------------------------------------+-------+-------+-------
      spatial function for dimension          |   0   |   1   |   2
     -----------------------------------------+-------+-------+-------
      origin of spatial function, x-component | x_[0] | x_[1] | x_[2]
     -----------------------------------------+-------+-------+-------
      origin of spatial function, y-component | y_[0] | y_[1] | y_[2]
     -----------------------------------------+-------+-------+-------
      origin of spatial function, x-component | z_[0] | z_[1] | z_[2]
     -----------------------------------------+-------+-------+-------

    */


    std::vector<double> x_; //! origin(s) of spatial function, x-component
    std::vector<double> y_; //! origin(s) of spatial function, y-component
    std::vector<double> z_; //! origin(s) of spatial function, z-component

    std::vector<Teuchos::RCP<DRT::PARSER::Parser<double> > > expr_; //! expression syntax tree(s)
  };


  /// special implementation for 3d Beltrami flow
  class BeltramiFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents()
      {
        return(4);
      };

  };


  /// special implementation for 2d Kim-Moin flow
  class KimMoinFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (this is a vector-valued functions)

    \return number of components (u,v,p)

    */
    virtual int NumberComponents()
      {
        return(3);
      };

  };


  /// special implementation for (randomly) disturbed 3d turbulent
  /// boundary-layer profile
  /// (currently fixed for low-Mach-number flow through a backward-facing step,
  ///  but may be easily manipulated to fit other profiles in other geometries)
  class TurbBouLayerFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents()
      {
        return(4);
      };

  };


  /// special implementation for Womersley blood flow
  class WomersleyFunction : public Function
  {
  public:

    /// ctor
    WomersleyFunction(bool locsys, int e, double radius, int mat, int curve);


    /// evaluate function at given position in space
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

  private:

    bool                   isinit_;
    bool                   locsys_;
    int                    locsysid_;
    double                 radius_;
    DRT::UTILS::TimeCurve& tc_;
    // exist after init phase:
    int                    mat_;
    int                    curve_;
    double                 viscosity_;
    double                 density_;

    // exist after init phase if locsys_==true
    // to do: move this stuff to separate class in locsys.H/.cpp
    // with functionality to transform spatial vectors from/to local system
    Condition*             locsyscond_;
    std::vector<double>    normal_;
    std::vector<double>    tangent1_;
    std::vector<double>    tangent2_;
    std::vector<double>    origin_;
  };


  /// special implementation for stationary 2d Jeffery-Hamel flow
  class JefferyHamelFlowFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents()
      {
        return(4);
      };

    /*!

    compute radial flow u as a function of the angle alpha

    \return radial velocity u

    \note this function is static such that we can get the radial velocity
          without creating the Function object

    */
    static double RadialVelocity(
        const double& alpha ///< angle between 0 and PI/4 (range is checked in debug mode)
        );

  };


  /// special implementation for a level set test function
  class ZalesaksDiskFunction : public Function
  {
  public:

    /// ctor
    ZalesaksDiskFunction();

    /// evaluate function at given position in space
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
  };

}
}


/*----------------------------------------------------------------------*/
//! Print function to be called from C
/*----------------------------------------------------------------------*/
extern "C"
void PrintFunctionDatHeader()
{
  DRT::UTILS::FunctionManager functionmanager;
  Teuchos::RCP<DRT::INPUT::Lines> lines = functionmanager.ValidFunctionLines();

  lines->Print(std::cout);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::INPUT::Lines> DRT::UTILS::FunctionManager::ValidFunctionLines()
{
  DRT::INPUT::LineDefinition linelin;
  linelin
    .AddNamedInt("FUNCT")
    .AddNamedDoubleVector("LINE_LIN",8)
    ;

  DRT::INPUT::LineDefinition linequad;
  linequad
    .AddNamedInt("FUNCT")
    .AddNamedDoubleVector("LINE_QUAD",6)
    ;

  DRT::INPUT::LineDefinition radiuslin;
  radiuslin
    .AddNamedInt("FUNCT")
    .AddNamedDoubleVector("RADIUS_LIN",8)
    ;

  DRT::INPUT::LineDefinition radiusquad;
  radiusquad
    .AddNamedInt("FUNCT")
    .AddNamedDoubleVector("RADIUS_QUAD",6)
    ;

  DRT::INPUT::LineDefinition beltrami;
  beltrami
    .AddNamedInt("FUNCT")
    .AddTag("BELTRAMI")
    ;

  DRT::INPUT::LineDefinition kimmoin;
  kimmoin
    .AddNamedInt("FUNCT")
    .AddTag("KIM-MOIN")
    ;

  DRT::INPUT::LineDefinition turbboulayer;
  turbboulayer
    .AddNamedInt("FUNCT")
    .AddTag("TURBBOULAYER")
    ;

  DRT::INPUT::LineDefinition jefferyhamel;
  jefferyhamel
    .AddNamedInt("FUNCT")
    .AddTag("JEFFERY-HAMEL")
    ;

  DRT::INPUT::LineDefinition womersley;
  womersley
    .AddNamedInt("FUNCT")
    .AddTag("WOMERSLEY")
    .AddNamedDouble("Radius")
    .AddNamedInt("MAT")
    .AddNamedInt("CURVE")
    ;

  DRT::INPUT::LineDefinition localwomersley;
  localwomersley
    .AddNamedInt("FUNCT")
    .AddTag("WOMERSLEY")
    .AddNamedInt("Local")
    .AddNamedDouble("Radius")
    .AddNamedInt("MAT")
    .AddNamedInt("CURVE")
    ;

  DRT::INPUT::LineDefinition cylinder3d;
  cylinder3d
    .AddNamedInt("FUNCT")
    .AddNamedDouble("CYLINDER_3D")
    ;

  DRT::INPUT::LineDefinition zalesaksdisk;
  zalesaksdisk
    .AddNamedInt("FUNCT")
    .AddTag("ZALESAKSDISK")
    ;

  DRT::INPUT::LineDefinition componentexpr;
  componentexpr
    .AddNamedInt("FUNCT")
    .AddNamedInt("COMPONENT")
    .AddNamedDoubleVector("EXPR",3)
    .AddNamedString("FUNCTION")
    ;

  DRT::INPUT::LineDefinition expr;
  expr
    .AddNamedInt("FUNCT")
    .AddNamedDoubleVector("EXPR",3)
    .AddNamedString("FUNCTION")
    ;

  Teuchos::RCP<DRT::INPUT::Lines> lines = Teuchos::rcp(new DRT::INPUT::Lines("FUNCT"));
  lines->Add(linelin);
  lines->Add(linequad);
  lines->Add(radiuslin);
  lines->Add(radiusquad);
  lines->Add(beltrami);
  lines->Add(kimmoin);
  lines->Add(turbboulayer);
  lines->Add(jefferyhamel);
  lines->Add(womersley);
  lines->Add(localwomersley);
  lines->Add(cylinder3d);
  lines->Add(zalesaksdisk);
  lines->Add(componentexpr);
  lines->Add(expr);
  return lines;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::FunctionManager::ReadInput(const DRT::INPUT::DatFileReader& reader)
{
  functions_.clear();

  Teuchos::RCP<DRT::INPUT::Lines> lines = ValidFunctionLines();

  // test for as many functions as there are
  for (int i = 1;; ++i)
  {
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition> > functions = lines->Read(reader, i);
    if (functions.size()==0)
      break;

    if (functions.size()==1)
    {
      Teuchos::RCP<DRT::INPUT::LineDefinition> function = functions[0];

      int id;
      function->ExtractInt("FUNCT",id);
      if (id!=i) dserror("expected FUNCT%d but got FUNCT%d", i, id);

      if (function->HaveNamed("LINE_LIN"))
      {
        std::vector<double> tmp;
        function->ExtractDoubleVector("LINE_LIN",tmp);

        double x1[3];
        double x2[3];
        x1[0] = tmp[0];
        x1[1] = tmp[1];
        x1[2] = tmp[2];
        //double val1 = tmp[3];
        x2[0] = tmp[4];
        x2[1] = tmp[5];
        x2[2] = tmp[6];
        //double val2 = tmp[7];

        /* calculate slope and offset */
        double b = tmp[3];
        double m = (tmp[7]-tmp[3]);
        double length = sqrt((tmp[4]-tmp[0])*(tmp[4]-tmp[0]) +
                             (tmp[5]-tmp[1])*(tmp[5]-tmp[1]) +
                             (tmp[6]-tmp[2])*(tmp[6]-tmp[2]));

        // Keep it simple. We use the expression based function class
        // that is able to handle straight lines.
        ostringstream expr;
        expr << "(" << b << ") + ((" << x2[0]-x1[0] << ")*x + (" << x2[1]-x1[1] << ")*y + (" << x2[2]-x1[2] << ")*z)/(" << length << ")/(" << length << ")*(" << m << ")";

        functions_.push_back(rcp(new ExprFunction(const_cast<char*>(expr.str().c_str()), x1[0], x1[1], x1[2])));
      }
      else if (function->HaveNamed("LINE_QUAD"))
      {
        std::vector<double> tmp;
        function->ExtractDoubleVector("LINE_QUAD",tmp);

        double x1[3];
        double x2[3];
        x1[0] = tmp[0];
        x1[1] = tmp[1];
        x1[2] = tmp[2];
        x2[0] = tmp[3];
        x2[1] = tmp[4];
        x2[2] = tmp[5];

        /* calculate length */
        double length = sqrt((tmp[3]-tmp[0])*(tmp[3]-tmp[0]) +
                             (tmp[4]-tmp[1])*(tmp[4]-tmp[1]) +
                             (tmp[5]-tmp[2])*(tmp[5]-tmp[2]));

        // Keep it simple.
        ostringstream expr;
        expr << "1.0 - 4 * (((" << x2[0]-x1[0] << ")*x + (" << x2[1]-x1[1] << ")*y + (" << x2[2]-x1[2] << ")*z)/(" << length << ")/(" << length << ") - 1.0/2.0)^2.0";
        functions_.push_back(rcp(new ExprFunction(const_cast<char*>(expr.str().c_str()), x1[0], x1[1], x1[2])));
      }
      else if (function->HaveNamed("RADIUS_LIN"))
      {
        std::vector<double> tmp;
        function->ExtractDoubleVector("RADIUS_LIN",tmp);

        double x1[3];
        double x2[3];
        x1[0] = tmp[0];
        x1[1] = tmp[1];
        x1[2] = tmp[2];
        //double val1 = tmp[3];
        x2[0] = tmp[4];
        x2[1] = tmp[5];
        x2[2] = tmp[6];
        //double val2 = tmp[7];

        /* calculate slope and offset */
        double b = tmp[3];
        double m = (tmp[7]-tmp[3]);
        double length = sqrt((tmp[4]-tmp[0])*(tmp[4]-tmp[0]) +
                             (tmp[5]-tmp[1])*(tmp[5]-tmp[1]) +
                             (tmp[6]-tmp[2])*(tmp[6]-tmp[2]));

        // Keep it simple.
        ostringstream expr;
        expr << "(" << b << ") + sqrt(x*x + y*y + z*z)/(" << length << ")*(" << m << ")";
        functions_.push_back(rcp(new ExprFunction(const_cast<char*>(expr.str().c_str()), x1[0], x1[1], x1[2])));
      }
      else if (function->HaveNamed("RADIUS_QUAD"))
      {
        std::vector<double> tmp;
        function->ExtractDoubleVector("RADIUS_QUAD",tmp);

        double x1[3];
        double x2[3];

        x1[0] = tmp[0];
        x1[1] = tmp[1];
        x1[2] = tmp[2];
        x2[0] = tmp[3];
        x2[1] = tmp[4];
        x2[2] = tmp[5];

        /* calculate length */
        double length = sqrt((tmp[3]-tmp[0])*(tmp[3]-tmp[0]) +
                             (tmp[4]-tmp[1])*(tmp[4]-tmp[1]) +
                             (tmp[5]-tmp[2])*(tmp[5]-tmp[2]));

        // Keep it simple.
        ostringstream expr;
        expr << "1.0 - (x*x + y*y + z*z)/(" << length << ")/(" << length << ")";
        functions_.push_back(rcp(new ExprFunction(const_cast<char*>(expr.str().c_str()), x1[0], x1[1], x1[2])));
      }
      else if (function->HaveNamed("BELTRAMI"))
      {
        functions_.push_back(rcp(new BeltramiFunction()));
      }
      else if (function->HaveNamed("KIM-MOIN"))
      {
        functions_.push_back(rcp(new KimMoinFunction()));
      }
      else if (function->HaveNamed("TURBBOULAYER"))
      {
        functions_.push_back(rcp(new TurbBouLayerFunction()));
      }
      else if (function->HaveNamed("JEFFERY-HAMEL"))
      {
        functions_.push_back(rcp(new JefferyHamelFlowFunction()));
      }
      else if (function->HaveNamed("WOMERSLEY"))
      {
        int e = -1;
        bool localcoordsystem = false;
        if (function->HaveNamed("Local"))
        {
          localcoordsystem = true;
          function->ExtractInt("Local",e);
        }

        double radius = -1.0;
        function->ExtractDouble("Radius", radius);
        int mat = -1;
        function->ExtractInt("MAT",mat);
        int curve = -1;
        function->ExtractInt("CURVE",curve);

        functions_.push_back(rcp(new WomersleyFunction(localcoordsystem,e-1,radius,mat,curve-1)));
      }
      else if (function->HaveNamed("CYLINDER_3D"))
      {
        double um;
        function->ExtractDouble("CYLINDER_3D", um);
        double h = 0.41;

        // Keep it simple.
        // This is particularly odd. Very special. Useless.
        ostringstream expr;
        expr << "16*(" << um << ")*y*z*((" << h << ")-y)*((" << h << ")-z) / ((" << h << ")^4)";
        functions_.push_back(rcp(new ExprFunction(const_cast<char*>(expr.str().c_str()), 0, 0, 0)));
      }
      else if (function->HaveNamed("ZALESAKSDISK"))
      {
        functions_.push_back(rcp(new ZalesaksDiskFunction()));
      }
      else if (function->HaveNamed("EXPR"))
      {
        Teuchos::RCP<ExprFunction> vecfunc = rcp(new ExprFunction());

        std::vector<double> origin;
        function->ExtractDoubleVector("EXPR",origin);
        std::string component;
        function->ExtractString("FUNCTION",component);

        vecfunc->AddExpr(component,origin[0],origin[1],origin[2]);
        functions_.push_back(vecfunc);
      }
      else
        dserror("unrecognized function");
    }

    else
    {
      Teuchos::RCP<ExprFunction> vecfunc = rcp(new ExprFunction());

      for (unsigned j=0; j<functions.size(); ++j)
      {
        int id;
        functions[j]->ExtractInt("FUNCT",id);
        if (id!=i) dserror("expected FUNCT%d but got FUNCT%d", i, id);

        if (not functions[j]->HaveNamed("COMPONENT"))
          dserror("component based expression function expected");

        int dim;
        functions[j]->ExtractInt("COMPONENT", dim);

        if (dim!=static_cast<int>(j))
        {
          dserror("For vector valued functions the components have to be\n"
                  "specified succesively, e.g. 0,1,..,ndof");
        }

        std::vector<double> origin;
        functions[j]->ExtractDoubleVector("EXPR",origin);
        std::string component;
        functions[j]->ExtractString("FUNCTION",component);

        vecfunc->AddExpr(component,origin[0],origin[1],origin[2]);
      }
      functions_.push_back(vecfunc);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::Function& DRT::UTILS::FunctionManager::Funct(int num)
{
  // ensure that desired function is available (prevents segmentation fault)
  if (functions_.size()< (unsigned int)(num+1) || num<0)
    dserror("function %d not available",num+1);

  return *(functions_[num]);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprFunction::ExprFunction()
{
  x_.clear();
  y_.clear();
  z_.clear();

  expr_.clear();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprFunction::ExprFunction(char* buf,
                                       double x,
                                       double y,
                                       double z)
{

  x_.push_back(x);
  y_.push_back(y);
  z_.push_back(z);

  expr_.push_back(rcp(new DRT::PARSER::Parser<double>(buf)));

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprFunction::~ExprFunction()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::ExprFunction::AddExpr(std::string buf,
                                       double x,
                                       double y,
                                       double z
  )
{
  expr_.push_back(rcp(new DRT::PARSER::Parser<double>(buf)));
  x_.push_back(x);
  y_.push_back(y);
  z_.push_back(z);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::ExprFunction::Evaluate(int index, const double* x, double t, DRT::Discretization* dis)
{
  // single expression for all components. Reset index to 0!
  if(expr_.size()==1)
  {
    index=0;
  }

  if(index>(int)expr_.size()-1 || index<0)
  {
    dserror("Tried to evaluate a function in a not available dimension.\nSpecify either one function or functions for all dimensions! \n(including one for the pressure)");
  }
  return expr_[index]->EvaluateFunct(x[0]-x_[index], x[1]-y_[index], x[2]-z_[index]);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::BeltramiFunction::Evaluate(int index, const double* xp, double t, DRT::Discretization* dis)
{
  double a = M_PI/4.0;
  double d = M_PI/2.0;

  switch (index)
  {
  case 0:
    return -a * ( exp(a*xp[0]) * sin(a*xp[1] + d*xp[2]) +
                  exp(a*xp[2]) * cos(a*xp[0] + d*xp[1]) );
  case 1:
    return -a * ( exp(a*xp[1]) * sin(a*xp[2] + d*xp[0]) +
                  exp(a*xp[0]) * cos(a*xp[1] + d*xp[2]) );
  case 2:
    return -a * ( exp(a*xp[2]) * sin(a*xp[0] + d*xp[1]) +
                  exp(a*xp[1]) * cos(a*xp[2] + d*xp[0]) );
  case 3:
    return -a*a/2 * ( exp(2*a*xp[0]) + exp(2*a*xp[1]) + exp(2*a*xp[2])
                      + 2* sin(a*xp[0]+d*xp[1]) * cos(a*xp[2]+d*xp[0]) * exp(a*(xp[1]+xp[2]))
                      + 2* sin(a*xp[1]+d*xp[2]) * cos(a*xp[0]+d*xp[1]) * exp(a*(xp[2]+xp[0]))
                      + 2* sin(a*xp[2]+d*xp[0]) * cos(a*xp[1]+d*xp[2]) * exp(a*(xp[0]+xp[1])));
  default:
    return 1.0;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::KimMoinFunction::Evaluate(int index, const double* xp, double t, DRT::Discretization* dis)
{
  double a = 2.0;

  switch (index)
  {
  case 0:
    return - cos(a*PI*xp[0]) * sin(a*PI*xp[1]);
  case 1:
    return + sin(a*PI*xp[0]) * cos(a*PI*xp[1]);
  case 2:
    return -1.0/4.0 * ( cos(2.0*a*PI*xp[0]) + cos(2.0*a*PI*xp[1]) );
  default:
    return 1.0;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::TurbBouLayerFunction::Evaluate(int index, const double* xp, double t, DRT::Discretization* dis)
{
  switch (index)
  {
    double randomnumber;
    double noise;

    case 0:
    {
      // set u_tau und nu (fixed for low-Mach-number flow through
      // backward-facing step, for the time being) and compute l_tau
      double utau = 0.106;
      double nu   = 0.00001527;
      double ltau = nu/utau;

      // compute y+ (upper wall of backward-facing step fixed, for the time being)
      double yplus = xp[1]/ltau;
      double myplus = (0.082-xp[1])/ltau;

      // generate noise via random number
      randomnumber = 2*((double)rand()-((double) RAND_MAX)/2.)/((double) RAND_MAX);
      noise = 0.1778 * randomnumber;

      // return velocity value in x-direction
      if      (yplus <= 5.0)
        return utau * yplus + noise;
      else if ((yplus > 5.0)   and (yplus <= 30.0))
        return utau * ( 5.0 * log(yplus) - 3.05 ) + noise;
      else if ((yplus > 30.0)  and (yplus <= 285.0))
      {
        double upre = utau * ( 2.5 * log(yplus) + 5.5 );
        if (upre > 2.063) return 2.063 + noise;
        else             return upre + noise;
      }
      else if ((myplus > 30.0) and (myplus <= 284.0))
      {
        double upre = utau * ( 2.5 * log(myplus) + 5.5 );
        if (upre > 2.063) return 2.063 + noise;
        else              return upre + noise;
      }
      else if ((myplus > 5.0)  and (myplus <= 30.0))
        return utau * ( 5.0 * log(myplus) - 3.05 ) + noise;
      else if (myplus <= 5.0)
        return utau * myplus + noise;
    }
    case 1:
    {
      // generate noise via random number
      randomnumber = 2*((double)rand()-((double) RAND_MAX)/2.)/((double) RAND_MAX);
      noise = 0.1778 * randomnumber;

      // return velocity value in y-direction
      return noise;
    }
    case 2:
    {
      // generate noise via random number
      randomnumber = 2*((double)rand()-((double) RAND_MAX)/2.)/((double) RAND_MAX);
      noise = 0.1778 * randomnumber;

      // return velocity value in z-direction
      return noise;
    }
    // nothing to be done for hydrodynamic pressure values
    case 3:
  default:
    return 0.0;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::JefferyHamelFlowFunction::RadialVelocity(
    const double& alpha
    )
{
#ifdef DEBUG
  if (alpha < 0.01 or alpha > (M_PI/4.0+0.01))
  {
    dserror("angle out of range! Must be between 0 and PI/4");
  }
#endif

  const double theta = alpha - (M_PI/8.0);

  // generated by Mathematica 6
  return -169.99995631379676
        + 51.72093368870616*    pow(theta, 2)
        + 1453.2480022627062*   pow(theta, 4)
        + 15469.249826734282*   pow(theta, 6)
        + 153515.23050166125*   pow(theta, 8)
        - 3.1393239563288596e6* pow(theta,10)
        + 1.3791074955303007e8* pow(theta,12)
        - 3.9054873809559045e9* pow(theta,14)
        + 8.064960114312076e10* pow(theta,16)
        - 1.2314175442399622e12*pow(theta,18)
        + 1.3952367263056582e13*pow(theta,20)
        - 1.1694967181678298e14*pow(theta,22)
        + 7.149819830278836e14* pow(theta,24)
        - 3.0970644442180215e15*pow(theta,26)
        + 9.001387683687223e15* pow(theta,28)
        - 1.5737446665792184e16*pow(theta,30)
        + 1.250536015803445e16* pow(theta,32);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::JefferyHamelFlowFunction::Evaluate(int index, const double* xp, double, DRT::Discretization*)
{

  const double x = xp[0];
  const double y = xp[1];

  const double alpha = atan(y/x);
  const double u_alpha = RadialVelocity(alpha);

  // you need to multiply with your viscosity outside of this function to
  // generate flow with higher Reynolds-numbers
  const double nu = 1;

  switch (index)
  {
  case 0:
    return  nu * (u_alpha/(x*x+y*y))*x;
  case 1:
    return  nu * (u_alpha/(x*x+y*y))*y;
  case 2:
    return 0.0;
  case 3:
    return 0.0;
  default:
    return 0.0;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::WomersleyFunction::WomersleyFunction(bool locsys, int e, double radius, int mat, int curve) :
Function(),
isinit_(false),
locsys_(locsys),
locsysid_(e),
radius_(radius),
tc_(DRT::Problem::Instance()->Curve(curve)),
mat_(mat),
curve_(curve),
viscosity_(-999.0e99),
density_(-999.0e99),
locsyscond_(NULL)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::UTILS::WomersleyFunction::Evaluate(int index, const double* xp, double t, DRT::Discretization* dis)
{
  if (!isinit_)
  {
    // get material parameters for fluid
    Teuchos::RCP<MAT::PAR::Material > mat = DRT::Problem::Instance()->Materials()->ById(mat_);
    if (mat->Type() != INPAR::MAT::m_fluid) dserror("Material %d is not a fluid",mat_);
    MAT::PAR::Parameter* params = mat->Parameter();
    MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
    if (!fparams) dserror("Material does not cast to Newtonian fluid");
    viscosity_ = fparams->viscosity_;
    density_ = fparams->density_;
    // get local coord system if any
    if (locsys_)
    {
      if (!dis) dserror("Have to pass a discretization to allow for local coord systems");
      vector<DRT::Condition*> locsys;
      dis->GetCondition("Locsys",locsys);
      if (!locsys.size()) dserror("No locsys conditions in discretization");
      for (int i=0; i<(int)locsys.size(); ++i)
        if (locsys[i]->Id() == locsysid_)
          locsyscond_ = locsys[i];
      if (!locsyscond_) dserror("Cannot find local coord system %d",locsysid_);
      const std::string* type  = locsyscond_->Get<std::string>("Type");
      if (*type != "FunctionEvaluation") dserror("Locsys is of wrong type %s",type->c_str());
      const vector<double>* n  = locsyscond_->Get<vector<double> >("normal");
      const vector<double>* t1 = locsyscond_->Get<vector<double> >("tangent");
      const vector<double>* o  = locsyscond_->Get<vector<double> >("origin");
      if (!n || !t1 || !o) dserror("Condition components missing");
      normal_   = *n;
      tangent1_ = *t1;
      origin_   = *o;
      tangent2_.resize(tangent1_.size());
      tangent2_[0] = normal_[1]*tangent1_[2]-normal_[2]*tangent1_[1];
      tangent2_[1] = normal_[2]*tangent1_[0]-normal_[0]*tangent1_[2];
      tangent2_[2] = normal_[0]*tangent1_[1]-normal_[1]*tangent1_[0];

      //printf("---------------------\n");
      //printf("n %f %f %f \nt1 %f %f %f \nt2 %f %f %f \no %f %f %f\n",
      //normal_[0],normal_[1],normal_[2],tangent1_[0],tangent1_[1],tangent1_[2],
      //tangent2_[0],tangent2_[1],tangent2_[2],origin_[0],origin_[1],origin_[2]);
    }
    isinit_ = true;
  }


  dserror("WomersleyFunction Evaluate not yet implemented");
  return -999.0;
}

/*----------------------------------------------------------------------*
 | constructor                                              henke 05/09 |
 *----------------------------------------------------------------------*/
DRT::UTILS::ZalesaksDiskFunction::ZalesaksDiskFunction() :
Function()
{
}

/*----------------------------------------------------------------------*
 | evaluation of level set test function "Zalesak's disk"   henke 05/09 |
 *----------------------------------------------------------------------*/
double DRT::UTILS::ZalesaksDiskFunction::Evaluate(int index, const double* xp, double t, DRT::Discretization* dis)
{
  //here calculation of distance (sign is already taken in consideration)
  double distance = 0;

  //inner part of slot
  if ((xp[0] <= 0.025) & (xp[0] >= -0.025) & (xp[1] >= 0.15) & (xp[1] <= 0.40))
  {
    distance = xp[1]-0.15;
    if ((xp[0]+0.025) < distance)
      distance = xp[0]+0.025;
    if ((0.025-xp[0]) < distance)
      distance = 0.025-xp[0];
  }

  //part directly under slot
  else if ((xp[0] <= 0.025) & (xp[0] >= -0.025) & (xp[1] >= 0.10) & (xp[1] < 0.15))
  {
    distance = xp[1]-0.15;
    if ((sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15) > distance)
      distance = sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15;
  }

  //part on the right hand side under the slot
  else if ((xp[0] > 0.025) & ((sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15)<0) & (xp[1] <= 0.15))
  {
    distance = sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15;
    if (-sqrt(DSQR(xp[0]-0.025)+DSQR(xp[1]-0.15)) > distance)
      distance = -sqrt(DSQR(xp[0]-0.025)+DSQR(xp[1]-0.15));
  }

  //part on the left hand side under the slot
  else if ((xp[0] < -0.025) & ((sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15) < 0) & (xp[1] <= 0.15))
  {
    distance = sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15;

    if (-sqrt(DSQR(xp[0]+0.025)+DSQR(xp[1]-0.15)) > distance)
      distance = -sqrt(DSQR(xp[0]+0.025)+DSQR(xp[1]-0.15));
  }

  //part on the right side of the slot
  else if ((xp[0]>0.025) & ((sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15)<0))
  {
    distance = sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15;
    if ((0.025-xp[0]) > distance)
      distance = 0.025-xp[0];
  }

  //part on the left side of the slot
  else if ((xp[0] < -0.025) & ((sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15)<0))
  {
    distance = sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15;
    if ((xp[0]+0.025) > distance)
      distance = xp[0]+0.025;
  }

  //trapezoidal part over the slot
  else if ((xp[1] > 0.40) & ((xp[1]-0.25-6*xp[0]) >= 0) & ((xp[1]-0.25+6*xp[0]) >= 0))
  {
    distance = sqrt(DSQR(xp[0]+0.025)+DSQR(xp[1]-0.4));
    if (sqrt(DSQR(xp[0]-0.025)+DSQR(xp[1]-0.40)) < distance)
      distance = sqrt(DSQR(xp[0]-0.025)+DSQR(xp[1]-0.40));
  }

  //rest of outer area
  else
  distance = sqrt(DSQR(xp[0])+DSQR(xp[1]-0.25))-0.15;

  return distance;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& DRT::UTILS::operator<<(std::ostream& out, const DRT::UTILS::Function& funct)
{
  out << "  Function:\n";
  return out;
}


#endif
