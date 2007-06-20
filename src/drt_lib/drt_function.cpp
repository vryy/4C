/*----------------------------------------------------------------------*/
/*!
\file drt_function.cpp

\brief Managing and evaluating of spatial functions

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>

#include "drt_function.H"

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*/
// the static instance
/*----------------------------------------------------------------------*/
DRT::FunctionManager DRT::FunctionManager::instance_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::FunctionManager::ReadInput()
{
  functions_.clear();

  // test for as many functions as there are
  for (int i = 1;; ++i)
  {
    ostringstream curve;
    curve << "--FUNCT" << i;
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
      frint("FUNCT",&id,&ierr);
      if (ierr!=1) dserror("cannot read FUNCT%d", i);
      if (id!=i) dserror("expected FUNCT%d but got FUNCT%d", i, id);

      /* read typ of funct */
      frchk("LINE_LIN",&ierr);
      if (ierr==1)
      {
        double tmp[8];
        frdouble_n("LINE_LIN",&(tmp[0]),8,&ierr);
        if (ierr!=1)
          dserror("failed to read function %d", i);

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

      frchk("LINE_QUAD",&ierr);
      if (ierr==1)
      {
        double tmp[6];
        frdouble_n("LINE_QUAD",&(tmp[0]),6,&ierr);

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

      frchk("RADIUS_LIN",&ierr);
      if (ierr==1)
      {
        double tmp[8];
        frdouble_n("RADIUS_LIN",&(tmp[0]),8,&ierr);

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

      frchk("RADIUS_QUAD",&ierr);
      if (ierr==1)
      {
        double tmp[6];
        frdouble_n("RADIUS_QUAD",&(tmp[0]),6,&ierr);

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

      frchk("BELTRAMI",&ierr);
      if (ierr==1)
      {
        functions_.push_back(rcp(new BeltramiFunction()));
      }

      frchk("KIM-MOIN",&ierr);
      if (ierr==1)
      {
        functions_.push_back(rcp(new KimMoinFunction()));
      }

      frchk("CYLINDER_3D",&ierr);
      if (ierr==1)
      {
        double um;
        frdouble_n("CYLINDER_3D",&um,1,&ierr);
        double h = 0.41;

        // Keep it simple.
        // This is particularly odd. Very special. Useless.
        ostringstream expr;
        expr << "16*(" << um << ")*y*z*((" << h << ")-y)*((" << h << ")-z) / ((" << h << ")^4)";
        functions_.push_back(rcp(new ExprFunction(const_cast<char*>(expr.str().c_str()), 0, 0, 0)));
      }

      frchk("EXPR",&ierr);
      if (ierr==1)
      {
        char expr[255];
        double x[3];

        /* read the position of the function's origin */
        frdouble_n("EXPR",x,3,&ierr);
        if (!ierr) dserror("failed to read coordinates");

        /* read the expression */
        frchar("FUNCTION", expr, &ierr);
        if (!ierr) dserror("failed to read expression string");
        functions_.push_back(rcp(new ExprFunction(expr, x[0], x[1], x[2])));
      }

#if 0
      frread();
      frchk("---",&ierr);
      if (ierr!=1)
        dserror("end of function definition expected");
#endif
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::ExprFunction::ExprFunction(char* buf, double x, double y, double z)
  : x_(x), y_(y), z_(z),
    expr_(pss_parse(buf))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::ExprFunction::~ExprFunction()
{
  pss_parse_cleanup(expr_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::ExprFunction::Evaluate(int index, const double* x)
{
  return pss_evaluate_funct(expr_, x[0]-x_, x[1]-y_, x[2]-z_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::BeltramiFunction::Evaluate(int index, const double* xp)
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
double DRT::KimMoinFunction::Evaluate(int index, const double* xp)
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
std::ostream& DRT::operator<<(std::ostream& out, const DRT::Function& funct)
{
  out << "  Function:\n";
  return out;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& DRT::operator<<(std::ostream& out, const DRT::FunctionManager& manager)
{
  out << "Function Manager:\n";
  for (unsigned i=0; i<manager.functions_.size(); ++i)
  {
    out << manager.functions_[i];
  }
  return out;
}

#endif
#endif
