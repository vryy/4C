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
#include "../drt_mat/newtonianfluid.H"


/*----------------------------------------------------------------------*/
// the static instance
/*----------------------------------------------------------------------*/
DRT::UTILS::FunctionManager DRT::UTILS::FunctionManager::instance_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::FunctionManager::ReadInput()
{
  functions_.clear();

  // test for as many functions as there are
  for (int i = 1;; ++i)
  {
    ostringstream curve;
    curve << "--FUNCT" << i;
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

      frchk("WOMERSLEY",&ierr);
      if (ierr==1)
      {
        int e = -1;
        bool localcoordsystem = false;
        frint("Local",&e,&ierr);
        if (ierr) localcoordsystem = true;
        double radius = -1.0;
        frdouble("Radius",&radius,&ierr);
        int mat = -1;
        frint("MAT",&mat,&ierr);
        if (!ierr) dserror("Word MAT missing in WOMERSLEY FUNCT");


        // input other stuff here....

        functions_.push_back(rcp(new WomersleyFunction(localcoordsystem,e-1,radius,mat)));
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

      frchk("ZALESAKSDISK",&ierr);
      if (ierr==1)
      {
        functions_.push_back(rcp(new ZalesaksDiskFunction(-1.0)));
      }

      frchk("EXPR",&ierr);
      if (ierr==1)
      {
        Teuchos::RCP<ExprFunction> vecfunc = rcp(new ExprFunction());

        for (int j=0;;++j)
        {
          char   component[255];
          double origin   [3];

          int    dim;
          frint("COMPONENT",&dim,&ierr);
          /* plausibility check */
          if (ierr == 1)
          {
            if(dim!=j)
            {
              dserror("For vector valued functions the components have to be \nspecified succesively, e.g. 0,1,..,ndof");
            }
          }

          /* read the position of the function's origin */
          frdouble_n("EXPR",origin,3,&ierr);
          if (!ierr) dserror("failed to read coordinates");

          /* read the expression */
          frchar("FUNCTION", component , &ierr);

          if (!ierr) dserror("failed to read expression string");

          (*vecfunc).AddExpr(component,origin[0],origin[1],origin[2]);

          frread();
          // stop if there is no content to this curve
          // no further curves are read
          frchk("---",&ierr);
          if (ierr==1)
          {
            break;
          }
        }

        functions_.push_back(vecfunc);

      }


#if 0
      frread();
      frchk("---",&ierr);
      if (ierr!=1)
        dserror("end of function definition expected");
#endif
    }
    else
    {
      // there is no such function, stop reading
      break;
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

  expr_.push_back(pss_parse(buf));

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprFunction::~ExprFunction()
{
  for(unsigned i=0;i!=expr_.size();++i)
  {
    pss_parse_cleanup(expr_[i]);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::ExprFunction::AddExpr(char* buf,
                                       double x,
                                       double y,
                                       double z
  )
{
  expr_.push_back(pss_parse(buf));
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
  return pss_evaluate_funct(expr_[index], x[0]-x_[index], x[1]-y_[index], x[2]-z_[index]);
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
DRT::UTILS::WomersleyFunction::WomersleyFunction(bool locsys, int e, double radius, int mat) :
Function(),
isinit_(false),
locsys_(locsys),
locsysid_(e),
radius_(radius),
mat_(mat),
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
      printf("---------------------\n");
      printf("n %f %f %f \nt1 %f %f %f \nt2 %f %f %f \no %f %f %f\n",
      normal_[0],normal_[1],normal_[2],tangent1_[0],tangent1_[1],tangent1_[2],
      tangent2_[0],tangent2_[1],tangent2_[2],origin_[0],origin_[1],origin_[2]);
    }
    isinit_ = true;
  }


  //dserror("WomersleyFunction Evaluate not yet implemented");
  return 0.0;
}

/*----------------------------------------------------------------------*
 | constructor                                              henke 05/09 |
 *----------------------------------------------------------------------*/
DRT::UTILS::ZalesaksDiskFunction::ZalesaksDiskFunction(double radius) :
Function(),
radius_(radius)
{
  /*
   * parameter "radius" is just an example for possible input parameters
   */
}

/*----------------------------------------------------------------------*
 | evaluation of level set test function "Zalesak's disk"   henke 05/09 |
 *----------------------------------------------------------------------*/
double DRT::UTILS::ZalesaksDiskFunction::Evaluate(int index, const double* xp, double t, DRT::Discretization* dis)
{
  dserror("Zalesak's disk function is not yet implemented");
  return 0.0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& DRT::UTILS::operator<<(std::ostream& out, const DRT::UTILS::Function& funct)
{
  out << "  Function:\n";
  return out;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::ostream& DRT::UTILS::operator<<(std::ostream& out, const DRT::UTILS::FunctionManager& manager)
{
  out << "Function Manager:\n";
  for (unsigned i=0; i<manager.functions_.size(); ++i)
  {
    out << manager.functions_[i];
  }
  return out;
}

#endif
