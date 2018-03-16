/*!-----------------------------------------------------------------------------------------------*
 \file fluid_functions.cpp

 \brief Managing and evaluating of spatial functions for fluid problems

  detailed description in header file

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "fluid_functions.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::BeltramiFunction::BeltramiFunction(double c1)
    : Function(),
      c1_(c1)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::BeltramiFunction::Evaluate(
    const int index,
    const double* xp,
    double t
    )
{
  double a = M_PI/4.0;
  double d = M_PI/2.0;

  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
  if (id==-1) dserror("Newtonian fluid material could not be found");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
  double dens = actmat->density_;
  double dynvisc = actmat->viscosity_;
  double kinvisc = dynvisc / dens;
  double tempfac = exp(-c1_*kinvisc*d*d*t);

  switch (index)
  {
  case 0:
    return -a * ( exp(a*xp[0]) * sin(a*xp[1] + d*xp[2]) +
                  exp(a*xp[2]) * cos(a*xp[0] + d*xp[1]) ) * tempfac;
  case 1:
    return -a * ( exp(a*xp[1]) * sin(a*xp[2] + d*xp[0]) +
                  exp(a*xp[0]) * cos(a*xp[1] + d*xp[2]) ) * tempfac;
  case 2:
    return -a * ( exp(a*xp[2]) * sin(a*xp[0] + d*xp[1]) +
                  exp(a*xp[1]) * cos(a*xp[2] + d*xp[0]) ) * tempfac;
  case 3:
    return -a*a/2 * dens * ( exp(2*a*xp[0]) + exp(2*a*xp[1]) + exp(2*a*xp[2])
                      + 2* sin(a*xp[0]+d*xp[1]) * cos(a*xp[2]+d*xp[0]) * exp(a*(xp[1]+xp[2]))
                      + 2* sin(a*xp[1]+d*xp[2]) * cos(a*xp[0]+d*xp[1]) * exp(a*(xp[2]+xp[0]))
                      + 2* sin(a*xp[2]+d*xp[0]) * cos(a*xp[1]+d*xp[2]) * exp(a*(xp[0]+xp[1]))) * tempfac;
  default:
    return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::BeltramiFunction::EvaluateTimeDerivative(const int      index,
                                                                  const double*  xp,
                                                                  const double   t,
                                                                  const unsigned deg)
{
  double a = M_PI/4.0;
  double d = M_PI/2.0;

  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
  if (id==-1) dserror("Newtonian fluid material could not be found");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
  double dens = actmat->density_;
  double dynvisc = actmat->viscosity_;
  double kinvisc = dynvisc / dens;
  double der1tempfac = (-c1_*kinvisc*d*d) * exp(-c1_*kinvisc*d*d*t);
  double der2tempfac = (-c1_*kinvisc*d*d) * (-c1_*kinvisc*d*d) * exp(-c1_*kinvisc*d*d*t);

  // resulting vector holding
  std::vector<double> res(deg+1);

  // add the value at time t
  res[0] = Evaluate(index,xp,t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    switch (index)
    {
    case 0:
      res[1] = -a * ( exp(a*xp[0]) * sin(a*xp[1] + d*xp[2]) +
                      exp(a*xp[2]) * cos(a*xp[0] + d*xp[1]) ) * der1tempfac;
    case 1:
      res[1] = -a * ( exp(a*xp[1]) * sin(a*xp[2] + d*xp[0]) +
                      exp(a*xp[0]) * cos(a*xp[1] + d*xp[2]) ) * der1tempfac;
    case 2:
      res[1] = -a * ( exp(a*xp[2]) * sin(a*xp[0] + d*xp[1]) +
                      exp(a*xp[1]) * cos(a*xp[2] + d*xp[0]) ) * der1tempfac;
    case 3:
      res[1] = -a*a/2 * dens * ( exp(2*a*xp[0]) + exp(2*a*xp[1]) + exp(2*a*xp[2])
                          + 2* sin(a*xp[0]+d*xp[1]) * cos(a*xp[2]+d*xp[0]) * exp(a*(xp[1]+xp[2]))
                          + 2* sin(a*xp[1]+d*xp[2]) * cos(a*xp[0]+d*xp[1]) * exp(a*(xp[2]+xp[0]))
                          + 2* sin(a*xp[2]+d*xp[0]) * cos(a*xp[1]+d*xp[2]) * exp(a*(xp[0]+xp[1]))) * der1tempfac;
    default:
      res[1] = 1.0;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    switch (index)
    {
    case 0:
      res[2] = -a * ( exp(a*xp[0]) * sin(a*xp[1] + d*xp[2]) +
                      exp(a*xp[2]) * cos(a*xp[0] + d*xp[1]) ) * der2tempfac;
    case 1:
      res[2] = -a * ( exp(a*xp[1]) * sin(a*xp[2] + d*xp[0]) +
                      exp(a*xp[0]) * cos(a*xp[1] + d*xp[2]) ) * der2tempfac;
    case 2:
      res[2] = -a * ( exp(a*xp[2]) * sin(a*xp[0] + d*xp[1]) +
                      exp(a*xp[1]) * cos(a*xp[2] + d*xp[0]) ) * der2tempfac;
    case 3:
      res[2] = -a*a/2 * dens * ( exp(2*a*xp[0]) + exp(2*a*xp[1]) + exp(2*a*xp[2])
                          + 2* sin(a*xp[0]+d*xp[1]) * cos(a*xp[2]+d*xp[0]) * exp(a*(xp[1]+xp[2]))
                          + 2* sin(a*xp[1]+d*xp[2]) * cos(a*xp[0]+d*xp[1]) * exp(a*(xp[2]+xp[0]))
                          + 2* sin(a*xp[2]+d*xp[0]) * cos(a*xp[1]+d*xp[2]) * exp(a*(xp[0]+xp[1]))) * der2tempfac;
    default:
      res[2] = 1.0;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::KimMoinFunction::Evaluate(const int index, const double* xp,
    double t)
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
std::vector<double> FLD::KimMoinFunction::EvaluateTimeDerivative(const int      index,
                                                                 const double*  xp,
                                                                 const double   t,
                                                                 const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg+1);

  // add the value at time t
  res[0] = Evaluate(index,xp,t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    switch (index)
    {
    case 0:
      res[1] = 0;
    case 1:
      res[1] = 0;
    case 2:
      res[1] = 0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    switch (index)
    {
    case 0:
      res[2] = 0;
    case 1:
      res[2] = 0;
    case 2:
      res[2] = 0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

return res;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::BochevUPFunction::Evaluate(const int index, const double* xp,
    double t)
{
  switch (index)
  {
  case 0:
    return sin(PI*xp[0]-0.7)*sin(PI*xp[1]+0.2);
  case 1:
    return cos(PI*xp[0]-0.7)*cos(PI*xp[1]+0.2);
  case 2:
    return sin(xp[0])*cos(xp[1])+(cos(1.0)-1.0)*sin(1.0);
  default:
    return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::BochevRHSFunction::Evaluate(const int index, const double* xp,
    double t)
{
  switch (index)
  {
  case 0:
    return cos(xp[0])*cos(xp[1])+PI*sin(PI*xp[0]-0.7)*cos(PI*xp[0]-0.7)+2*PI*PI*sin(PI*xp[0]-0.7)*sin(PI*xp[1]+0.2);
  case 1:
    return -sin(xp[0])*sin(xp[1])-PI*sin(PI*xp[1]+0.2)*cos(PI*xp[1]+0.2)+2*PI*PI*cos(PI*xp[0]-0.7)*cos(PI*xp[1]+0.2);
  default:
    return 1.0;
  }
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::BeltramiUP::BeltramiUP(int mat_id) :
Function(),
density_(-999.0e99),
kinviscosity_(-999.0e99)
{

  // get material parameters for fluid
  Teuchos::RCP<MAT::PAR::Material > mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
  if (mat->Type() != INPAR::MAT::m_fluid)
    dserror("Material %d is not a fluid",mat_id);
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get density
  density_ = fparams->density_;

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / density_;

}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::BeltramiUP::BeltramiUP( Teuchos::RCP<MAT::Material> & mat) :
Function(),
density_(-999.0e99),
kinviscosity_(-999.0e99)
{

  if (mat->MaterialType() != INPAR::MAT::m_fluid)
    dserror("Material is not a fluid");
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get density
  density_ = fparams->density_;

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / fparams->density_;

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::BeltramiUP::Evaluate(const int index, const double* xp,
    double t)
{

  double x = xp[0];
  double y = xp[1];
  double z = xp[2];

  double a = PI/4.0;
  double b = PI/4.0;
  double c= a*a + b*b + a*b;

  double K1 = exp(a*(x-z) + b*(y-z)); // =K4
  double K2 = exp(a*(z-y) + b*(x-y)); // =K5
  double K3 = exp(a*(y-x) + b*(z-x)); // =K6


  switch (index)
  {
  case 0:
    return b*K1-a*K2;
  case 1:
    return b*K3-a*K1;
  case 2:
    return b*K2-a*K3;
  case 3:
    return  c*( 1.0/K3 + 1.0/K2 + 1.0/K1 )* density_;
  default:
    dserror("wrong index %d", index);
    break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::BeltramiUP::EvaluateTimeDerivative(const int      index,
                                                            const double*  xp,
                                                            const double   t,
                                                            const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg+1);

  // add the value at time t
  res[0] = Evaluate(index,xp,t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    switch (index)
    {
    case 0:
      res[1] = 0;
    case 1:
      res[1] = 0;
    case 2:
      res[1] = 0;
    case 3:
      res[1] = 0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    switch (index)
    {
    case 0:
      res[2] = 0;
    case 1:
      res[2] = 0;
    case 2:
      res[2] = 0;
    case 3:
      res[2] = 0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

return res;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::BeltramiGradU::BeltramiGradU(int mat_id ) :
Function(),
kinviscosity_(-999.0e99)
{

  // get material parameters for fluid
  Teuchos::RCP<MAT::PAR::Material > mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
  if (mat->Type() != INPAR::MAT::m_fluid)
    dserror("Material %d is not a fluid",mat_id);
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / fparams->density_;

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::BeltramiGradU::BeltramiGradU(Teuchos::RCP<MAT::Material> & mat ) :
Function(),
kinviscosity_(-999.0e99)
{

  if (mat->MaterialType() != INPAR::MAT::m_fluid)
    dserror("Material is not a fluid");
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / fparams->density_;

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::BeltramiGradU::Evaluate(const int index, const double* xp,
    double t)
{

  double x = xp[0];
  double y = xp[1];
  double z = xp[2];

  double a = PI/4.0;
  double b = PI/4.0;

  double K1 = exp(a*(x-z) + b*(y-z)); // =K4
  double K2 = exp(a*(z-y) + b*(x-y)); // =K5
  double K3 = exp(a*(y-x) + b*(z-x)); // =K6


  switch (index)
  {
  case 0: // u,x
    return a*b*(K1-K2);
  case 1: // u,y
    return b*b*K1+a*(a+b)*K2;
  case 2: // u,z
    return -b*(a+b)*K1-a*a*K2;
  case 3: // v,x
    return  -b*(a+b)*K3-a*a*K1;
  case 4: // v,y
    return a*b*K3-a*b*K1;
  case 5: // v,z
    return b*b*K3+a*(a+b)*K1;
  case 6: // w,x
    return b*b*K2+a*(a+b)*K3;
  case 7: // w,y
    return -b*(a+b)*K2-a*a*K3;
  case 8: // w,z
    return a*b*(K2-K3);
  default:
    dserror("wrong index %d", index);
    break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::BeltramiGradU::EvaluateTimeDerivative(const int      index,
                                                            const double*  xp,
                                                            const double   t,
                                                            const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg+1);

  // add the value at time t
  res[0] = Evaluate(index,xp,t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    switch (index)
    {
    case 0:
      res[1] = 0;
    case 1:
      res[1] = 0;
    case 2:
      res[1] = 0;
    case 3:
      res[1] = 0;
    case 4:
      res[1] = 0;
    case 5:
      res[1] = 0;
    case 6:
      res[1] = 0;
    case 7:
      res[1] = 0;
    case 8:
      res[1] = 0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    switch (index)
    {
    case 0:
      res[2] = 0;
    case 1:
      res[2] = 0;
    case 2:
      res[2] = 0;
    case 3:
      res[2] = 0;
    case 4:
      res[2] = 0;
    case 5:
      res[2] = 0;
    case 6:
      res[2] = 0;
    case 7:
      res[2] = 0;
    case 8:
      res[2] = 0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

return res;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::BeltramiRHS::BeltramiRHS(int mat_id, bool is_stokes) :
Function(),
kinviscosity_(-999.0e99),
is_stokes_(is_stokes)
{

  // get material parameters for fluid
  Teuchos::RCP<MAT::PAR::Material > mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
  if (mat->Type() != INPAR::MAT::m_fluid)
    dserror("Material %d is not a fluid",mat_id);
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / fparams->density_;

}


double FLD::BeltramiRHS::Evaluate(const int index, const double* xp,
    double t)
{

  double x = xp[0];
  double y = xp[1];
  double z = xp[2];

  double a = PI/4.0;
  double b = PI/4.0;
  double c= a*a + b*b + a*b;

  double K1 = exp(a*(x-z) + b*(y-z)); // =K4
  double K2 = exp(a*(z-y) + b*(x-y)); // =K5
  double K3 = exp(a*(y-x) + b*(z-x)); // =K6

  double t1 = b*K1-a*K2;
  double t2 = b*K3-a*K1;
  double t3 = b*K2-a*K3;

  double conv_x = 0.0;
  double conv_y = 0.0;
  double conv_z = 0.0;

  if(!is_stokes_)
  {
    conv_x = t1 * (     a*b*K1 -     a*b*K2) + t2 * (     b*b*K1 + a*(a+b)*K2) + t3 * (-b*(a+b)*K1 -     a*a*K2);
    conv_y = t1 * (-b*(a+b)*K3 -     a*a*K1) + t2 * (     a*b*K3 -     a*b*K1) + t3 * (     b*b*K3 + a*(a+b)*K1);
    conv_z = t1 * (     b*b*K2 + a*(a+b)*K3) + t2 * (-b*(a+b)*K2 -     a*a*K3) + t3 * (     a*b*K2 -     a*b*K3);
  }

  switch (index)
  {
  case 0:
    return c*( (a+b)/K3 - b/K2 - a/K1 - 2.*kinviscosity_*(b*K1-a*K2) ) + conv_x;
  case 1:
    return c*( -a/K3 + (a+b)/K2 - b/K1 - 2.*kinviscosity_*(b*K3-a*K1) ) + conv_y;
  case 2:
    return c*( -b/K3 - a/K2 + (a+b)/K1 - 2.*kinviscosity_*(b*K2-a*K3) ) + conv_z;
  default:
    dserror("wrong index %d", index);
    break;
  }

  return 1.0;

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::BeltramiRHS::EvaluateTimeDerivative(const int      index,
                                                            const double*  xp,
                                                            const double   t,
                                                            const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg+1);

  // add the value at time t
  res[0] = Evaluate(index,xp,t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    switch (index)
    {
    case 0:
      res[1] = 0;
    case 1:
      res[1] = 0;
    case 2:
      res[1] = 0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    switch (index)
    {
    case 0:
      res[2] = 0;
    case 1:
      res[2] = 0;
    case 2:
      res[2] = 0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

return res;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::KimMoinUP::KimMoinUP(int mat_id, bool is_stationary) :
Function(),
density_(-999.0e99),
kinviscosity_(-999.0e99),
is_stationary_(is_stationary)
{

  // get material parameters for fluid
  Teuchos::RCP<MAT::PAR::Material > mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
  if (mat->Type() != INPAR::MAT::m_fluid)
    dserror("Material %d is not a fluid",mat_id);
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get density
  density_ = fparams->density_;

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / density_;

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::KimMoinUP::KimMoinUP( Teuchos::RCP<MAT::Material> & mat,
    bool is_stationary) :
Function(),
density_(-999.0e99),
kinviscosity_(-999.0e99),
is_stationary_(is_stationary)
{

  if (mat->MaterialType() != INPAR::MAT::m_fluid)
    dserror("Material is not a fluid");
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get density
  density_ = fparams->density_;

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / fparams->density_;

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::KimMoinUP::Evaluate(const int index, const double* xp,
    double t)
{

  double x = xp[0];
  double y = xp[1];
  //double z = xp[2];

  double a = 2.0;

  double a_pi_x = a*PI*x;
  double a_pi_y = a*PI*y;

  // time factors
  double gu = 1.0;
  double gp = 1.0;

  if(!is_stationary_)
  {
    gu = exp(-2.0*a*a*PI*PI*t*kinviscosity_);
    gp = exp(-4.0*a*a*PI*PI*t*kinviscosity_);
  }

  switch (index)
  {
  case 0:
    return -cos(a_pi_x)*sin(a_pi_y) * gu;
  case 1:
    return sin(a_pi_x)*cos(a_pi_y) * gu;
  case 2:
    return 0.0;
  case 3:
    return  -1./4. * ( cos(2.0*a_pi_x) + cos(2.0*a_pi_y) ) * gp * density_;
  default:
    dserror("wrong index %d", index);
    break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::KimMoinUP::EvaluateTimeDerivative(const int      index,
                                                           const double*  xp,
                                                           const double   t,
                                                           const unsigned deg)
{
  double x = xp[0];
  double y = xp[1];
  //double z = xp[2];

  double a = 2.0;

  double a_pi_x = a*PI*x;
  double a_pi_y = a*PI*y;

  // resulting vector holding
  std::vector<double> res(deg+1);

  // add the value at time t
  res[0] = Evaluate(index,xp,t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // time factors
    double gu = 1.0;
    double gp = 1.0;

    if(!is_stationary_)
    {
      gu = (-2.0*a*a*PI*PI*kinviscosity_)*exp(-2.0*a*a*PI*PI*t*kinviscosity_);
      gp = (-4.0*a*a*PI*PI*kinviscosity_)*exp(-4.0*a*a*PI*PI*t*kinviscosity_);
    }

    switch (index)
    {
    case 0:
      res[1] = -cos(a_pi_x)*sin(a_pi_y) * gu;
    case 1:
      res[1] = sin(a_pi_x)*cos(a_pi_y) * gu;
    case 2:
      res[1] = 0.0;
    case 3:
      res[1] = -1./4. * ( cos(2.0*a_pi_x) + cos(2.0*a_pi_y) ) * gp * density_;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    // time factors
    double gu = 1.0;
    double gp = 1.0;

    if(!is_stationary_)
    {
      gu = (-2.0*a*a*PI*PI*kinviscosity_)*(-2.0*a*a*PI*PI*kinviscosity_)*exp(-2.0*a*a*PI*PI*t*kinviscosity_);
      gp = (-4.0*a*a*PI*PI*kinviscosity_)*(-4.0*a*a*PI*PI*kinviscosity_)*exp(-4.0*a*a*PI*PI*t*kinviscosity_);
    }

    switch (index)
    {
    case 0:
      res[2] = -cos(a_pi_x)*sin(a_pi_y) * gu;
    case 1:
      res[2] =  sin(a_pi_x)*cos(a_pi_y) * gu;
    case 2:
      res[2] = 0.0;
    case 3:
      res[2] = -1./4. * ( cos(2.0*a_pi_x) + cos(2.0*a_pi_y) ) * gp * density_;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

return res;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::KimMoinGradU::KimMoinGradU(int mat_id, bool is_stationary ) :
Function(),
kinviscosity_(-999.0e99),
is_stationary_(is_stationary)
{

  // get material parameters for fluid
  Teuchos::RCP<MAT::PAR::Material > mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
  if (mat->Type() != INPAR::MAT::m_fluid)
    dserror("Material %d is not a fluid",mat_id);
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / fparams->density_;

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::KimMoinGradU::KimMoinGradU(Teuchos::RCP<MAT::Material> & mat, bool is_stationary ) :
Function(),
kinviscosity_(-999.0e99),
is_stationary_(is_stationary)
{

  if (mat->MaterialType() != INPAR::MAT::m_fluid)
    dserror("Material is not a fluid");
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / fparams->density_;

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::KimMoinGradU::Evaluate(const int index, const double* xp,
    double t)
{
  //double visc = 1.0;

  double x = xp[0];
  double y = xp[1];
  //double z = xp[2];

  double a = 2.0;

  double a_pi_x = a*PI*x;
  double a_pi_y = a*PI*y;

  // time factors
  double gu = 1.0;

  if(!is_stationary_)
  {
    gu = exp(-2.0*a*a*PI*PI*t*kinviscosity_);
  }


  switch (index)
  {
  case 0: // u,x
    return  sin(a_pi_x)*sin(a_pi_y)*a*PI * gu;
  case 1: // u,y
    return -cos(a_pi_x)*cos(a_pi_y)*a*PI * gu;
  case 2: // u,z
    return 0.0;
  case 3: // v,x
    return  cos(a_pi_x)*cos(a_pi_y)*a*PI * gu;
  case 4: // v,y
    return -sin(a_pi_x)*sin(a_pi_y)*a*PI * gu;
  case 5: // v,z
    return 0.0;
  case 6: // w,x
    return 0.0;
  case 7: // w,y
    return 0.0;
  case 8: // w,z
    return 0.0;
  default:
    dserror("wrong index %d", index);
    break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::KimMoinGradU::EvaluateTimeDerivative(const int      index,
                                                              const double*  xp,
                                                              const double   t,
                                                              const unsigned deg)
{
  //double visc = 1.0;

  double x = xp[0];
  double y = xp[1];
  //double z = xp[2];

  double a = 2.0;

  double a_pi_x = a*PI*x;
  double a_pi_y = a*PI*y;

  // resulting vector holding
  std::vector<double> res(deg+1);

  // add the value at time t
  res[0] = Evaluate(index,xp,t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // time factors
    double gu = 1.0;

    if(!is_stationary_)
    {
      gu = (-2.0*a*a*PI*PI*kinviscosity_)*exp(-2.0*a*a*PI*PI*t*kinviscosity_);
    }

    switch (index)
    {
    case 0: // u,x
      res[1] =  sin(a_pi_x)*sin(a_pi_y)*a*PI * gu;
    case 1: // u,y
      res[1] = -cos(a_pi_x)*cos(a_pi_y)*a*PI * gu;
    case 2: // u,z
      res[1] = 0.0;
    case 3: // v,x
      res[1] =  cos(a_pi_x)*cos(a_pi_y)*a*PI * gu;
    case 4: // v,y
      res[1] = -sin(a_pi_x)*sin(a_pi_y)*a*PI * gu;
    case 5: // v,z
      res[1] = 0.0;
    case 6: // w,x
      res[1] = 0.0;
    case 7: // w,y
      res[1] = 0.0;
    case 8: // w,z
      res[1] = 0.0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    // time factors
    double gu = 1.0;

    if(!is_stationary_)
    {
      gu = (-2.0*a*a*PI*PI*kinviscosity_)*(-2.0*a*a*PI*PI*kinviscosity_)*exp(-2.0*a*a*PI*PI*t*kinviscosity_);
    }

    switch (index)
    {
    case 0: // u,x
      res[2] =  sin(a_pi_x)*sin(a_pi_y)*a*PI * gu;
    case 1: // u,y
      res[2] = -cos(a_pi_x)*cos(a_pi_y)*a*PI * gu;
    case 2: // u,z
      res[2] = 0.0;
    case 3: // v,x
      res[2] =  cos(a_pi_x)*cos(a_pi_y)*a*PI * gu;
    case 4: // v,y
      res[2] = -sin(a_pi_x)*sin(a_pi_y)*a*PI * gu;
    case 5: // v,z
      res[2] = 0.0;
    case 6: // w,x
      res[2] = 0.0;
    case 7: // w,y
      res[2] = 0.0;
    case 8: // w,z
      res[2] = 0.0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

return res;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::KimMoinRHS::KimMoinRHS(int mat_id, bool is_stationary,
    bool is_stokes) :
Function(),
kinviscosity_(-999.0e99),
is_stationary_(is_stationary),
is_stokes_(is_stokes)
{

  // get material parameters for fluid
  Teuchos::RCP<MAT::PAR::Material > mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
  if (mat->Type() != INPAR::MAT::m_fluid)
    dserror("Material %d is not a fluid",mat_id);
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / fparams->density_;

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::KimMoinRHS::Evaluate(const int index, const double* xp,
    double t)
{
  double x = xp[0];
  double y = xp[1];
  //double z = xp[2];

  double a = 2.0;

  // time factors
  double gu = 1.0;
  double gp = 1.0;

  if(!is_stationary_)
  {
    gu = exp(-2.0*a*a*PI*PI*t*kinviscosity_);
    gp = exp(-4.0*a*a*PI*PI*t*kinviscosity_);
  }


  double a_pi_x = a*PI*x;
  double a_pi_y = a*PI*y;

  double conv_x = 0.0;
  double conv_y = 0.0;
  //double conv_z = 0.0;

  if(!is_stokes_)
  {
    conv_x = -a*PI*sin(a_pi_x)*cos(a_pi_x);
    conv_y = -a*PI*sin(a_pi_y)*cos(a_pi_y);
  }

  double visc_x = 0.0;
  double visc_y = 0.0;

  if(is_stationary_)
  {
    visc_x = kinviscosity_*2.*a*a*PI*PI* (-cos(a_pi_x)*sin(a_pi_y)); // * (gu = 1) for stationary
    visc_y = kinviscosity_*2.*a*a*PI*PI* ( sin(a_pi_x)*cos(a_pi_y)); // * (gu = 1) for stationary
  }

  // in case of instationary: du/dt - \nu \laplacian(u) = 0

  switch (index)
  {
  case 0:
    return 0.5*a*PI*sin(2.*a_pi_x)*gp + visc_x + conv_x*gu*gu;
  case 1:
    return 0.5*a*PI*sin(2.*a_pi_y)*gp + visc_y + conv_y*gu*gu;
  case 2:
    return 0.0;
  default:
    dserror("wrong index %d", index);
    break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::KimMoinRHS::EvaluateTimeDerivative(const int      index,
                                                            const double*  xp,
                                                            const double   t,
                                                            const unsigned deg)
{
  double x = xp[0];
  double y = xp[1];
  //double z = xp[2];

  double a = 2.0;

  // resulting vector holding
  std::vector<double> res(deg+1);

  // add the value at time t
  res[0] = Evaluate(index,xp,t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // time factors
    double gu = 1.0;
    double gp = 1.0;

    if(!is_stationary_)
    {
      gu = (-4.0*a*a*PI*PI*kinviscosity_)*exp(-4.0*a*a*PI*PI*t*kinviscosity_);
      gp = (-4.0*a*a*PI*PI*kinviscosity_)*exp(-4.0*a*a*PI*PI*t*kinviscosity_);
    }

    double a_pi_x = a*PI*x;
    double a_pi_y = a*PI*y;

    double conv_x = 0.0;
    double conv_y = 0.0;
    //double conv_z = 0.0;

    if(!is_stokes_)
    {
      conv_x = -a*PI*sin(a_pi_x)*cos(a_pi_x);
      conv_y = -a*PI*sin(a_pi_y)*cos(a_pi_y);
    }

    // in case of instationary: du/dt - \nu \laplacian(u) = 0

    switch (index)
    {
    case 0:
      res[1] = 0.5*a*PI*sin(2.*a_pi_x)*gp + conv_x*gu;
    case 1:
      res[1] = 0.5*a*PI*sin(2.*a_pi_y)*gp + conv_y*gu;
    case 2:
      res[1] = 0.0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    // time factors
    double gu = 1.0;
    double gp = 1.0;

    if(!is_stationary_)
    {
      gu = (-4.0*a*a*PI*PI*kinviscosity_)*(-4.0*a*a*PI*PI*kinviscosity_)*exp(-4.0*a*a*PI*PI*t*kinviscosity_);
      gp = (-4.0*a*a*PI*PI*kinviscosity_)*(-4.0*a*a*PI*PI*kinviscosity_)*exp(-4.0*a*a*PI*PI*t*kinviscosity_);
    }

    double a_pi_x = a*PI*x;
    double a_pi_y = a*PI*y;

    double conv_x = 0.0;
    double conv_y = 0.0;
    //double conv_z = 0.0;

    if(!is_stokes_)
    {
      conv_x = -a*PI*sin(a_pi_x)*cos(a_pi_x);
      conv_y = -a*PI*sin(a_pi_y)*cos(a_pi_y);
    }

    // in case of instationary: du/dt - \nu \laplacian(u) = 0

    switch (index)
    {
    case 0:
      res[2] = 0.5*a*PI*sin(2.*a_pi_x)*gp + conv_x*gu;
    case 1:
      res[2] = 0.5*a*PI*sin(2.*a_pi_y)*gp + conv_y*gu;
    case 2:
      res[2] = 0.0;
    default:
      dserror("wrong index %d", index);
      break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

return res;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::KimMoinStress::KimMoinStress(int mat_id, bool is_stationary, double amplitude) :
Function(),
kinviscosity_(-999.0e99),
is_stationary_(is_stationary),
amplitude_(amplitude)
{

  if (amplitude_ != 1.0)
    dserror("At the moment other implementation of Kimmoin Functions do not include the amplitude functionality!");

  // get material parameters for fluid
  Teuchos::RCP<MAT::PAR::Material > mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
  if (mat->Type() != INPAR::MAT::m_fluid)
    dserror("Material %d is not a fluid",mat_id);
  MAT::PAR::Parameter* params = mat->Parameter();
  MAT::PAR::NewtonianFluid* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
  if (!fparams)
    dserror("Material does not cast to Newtonian fluid");

  // get density
  density_ = fparams->density_;

  // get kinematic viscosity
  kinviscosity_ = fparams->viscosity_ / density_;


}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::KimMoinStress::Evaluate(int index, const double* xp, double t)
{
  //double visc = 1.0;

  double x = xp[0];
  double y = xp[1];
  //double z = xp[2];

  const double a = 2;// in general does not need to be equal ... PARAM_A;
  const double b = 2;// in general does not need to be equal ... PARAM_B;

  double a_pi_x = a*PI*x;
  double a_pi_y = a*PI*y;

  // time factors
  double gu = 1.0;
  double gp = 1.0;

  if(!is_stationary_)
  {
    gu = exp(-2.0*b*b*PI*PI*t*kinviscosity_);
    gp = exp(-4.0*b*b*PI*PI*t*kinviscosity_);
  }

  double fac = 2*kinviscosity_*density_*gu*PI*a*sin(a_pi_x)*sin(a_pi_y)*amplitude_;
  double p = -1./4. * ( cos(2.0*a_pi_x) + cos(2.0*a_pi_y) ) * gp * density_;


  switch (index)
  {
  case 0: // sigma_xx
    return  (fac - p);
  case 1: // sigma_yy
    return (-fac - p);
  case 2: // sigma_zz
    return (-p);
  case 3: // sigma_xy
    return  0.0;
  case 4: // sigma_yz
    return 0.0;
  case 5: // sigma_zx
    return 0.0;
  default:
    dserror("wrong index %d", index);
    break;
  }

  return 1.0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::TurbBouLayerFunction::Evaluate(const int index, const double* xp,
    double t)
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
      randomnumber = DRT::Problem::Instance()->Random()->Uni();
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
      randomnumber = DRT::Problem::Instance()->Random()->Uni();
      noise = 0.1778 * randomnumber;

      // return velocity value in y-direction
      return noise;
    }
    case 2:
    {
      // generate noise via random number
      randomnumber = DRT::Problem::Instance()->Random()->Uni();
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
double FLD::TurbBouLayerFunctionBFS::Evaluate(const int index,
    const double* xp, double t)
{
  double randomnumber = 0.0;
  double noise = 0.0;
  // maximal velocity
  const double max = 1.899;
  // fluctuation 10% of maximal velocity
  const double fluc = 0.1 * max;
  // generate noise in turbulent boundary layer only
  const bool boulayer = true;
  // step height
  double h = 0.041;
  // boundary layer thickness
  double delta = 1.2 * h;

  switch (index)
  {
    case 0:
    {
      // set u_tau und nu and compute l_tau
      double utau = 0.093;
      double nu   = 0.000015268;
      double ltau = nu/utau;

      // compute y+
      double yplus = xp[1]/ltau;

      if (!boulayer)
      {
        // noise over complete inlet
        // generate noise via random number between -1 and 1
        randomnumber = DRT::Problem::Instance()->Random()->Uni();
        noise = fluc * randomnumber;
      }
      else
      {
        // noise only in boundary layer
        if (xp[1]<=delta)
        {
          // generate noise via random number between -1 and 1
          randomnumber = DRT::Problem::Instance()->Random()->Uni();
          noise = fluc * randomnumber;
        }
        else
        {
          noise = 0.0;
        }
      }

      // return velocity value in x-direction
      // viscous sublayer
      if (yplus <= 5.0)
        return utau * yplus + noise;
      // buffer layer
      else if ((yplus > 5.0) and (yplus <= 30.0))
        return utau * ( 4.5 * log(yplus) - 2.2 ) + noise;
      // log-layer
      else if (yplus > 30.0 and (yplus <= 150.0))
        return utau * ( 2.6 * log(yplus) + 4.3 ) + noise;
      // defect law
      else if (yplus > 150.0)
      {
        double upre = utau * ( 3.6 * log(xp[1]/delta) - 0.35 ) + max;
        if (upre > max) return max + noise;
        else            return upre + noise;
      }
    }
    case 1:
    {
      if (!boulayer)
      {
        // noise over complete inlet
        // generate noise via random number between -1 and 1
        randomnumber = DRT::Problem::Instance()->Random()->Uni();
        noise = fluc * randomnumber;
      }
      else
      {
        // noise only in boundary layer
        if (xp[1]<=delta)
        {
          // generate noise via random number between -1 and 1
          randomnumber = DRT::Problem::Instance()->Random()->Uni();
          noise = fluc * randomnumber;
        }
        else
        {
          noise = 0.0;
        }
      }

      // return velocity value in y-direction
      return noise;
    }
    case 2:
    {
      if (!boulayer)
      {
        // noise over complete inlet
        // generate noise via random number between -1 and 1
        randomnumber = DRT::Problem::Instance()->Random()->Uni();
        noise = fluc * randomnumber;
      }
      else
      {
        // noise only in boundary layer
        if (xp[1]<=delta)
        {
          // generate noise via random number between -1 and 1
         randomnumber = DRT::Problem::Instance()->Random()->Uni();
          noise = fluc * randomnumber;
        }
        else
        {
          noise = 0.0;
        }
      }

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
double FLD::TurbBouLayerFunctionORACLES::Evaluate(const int index,
    const double* xp, double t)
{
  double randomnumber;
  double noise;

  // fluctuation 10% of bulk mean velocity
  double fluc= 0.1;
  // maximal velocity
  double umax = 15.62;

  switch (index)
  {
    case 0:
    {
      // set u_tau und nu and compute l_tau
      const double utau = 0.7;
      const double nu   = 1.375E-5;
      double ltau = nu/utau;

      // transform global y-coordinate to local (upper and lower) channel y-coordinates
      double y=0.0;
      // remark: no tolerances needed, since there will always be no-slip boundary conditions
      // upper half upper inlet channel
      if (xp[1]>0.0202 and xp[1]<0.0354+1.0E-9)
        y = -(xp[1]-0.0354);
      // lower half upper inlet channel
      else if (xp[1]>0.005-1.0E-9 and xp[1]<=0.0202)
        y = xp[1]-0.005;
      // upper half lower inlet channel
      else if (xp[1]>=-0.0202 and xp[1]<-0.005+1.0E-9)
        y = -(xp[1]+0.005);
      // lower half lower inlet channel
      else if (xp[1]>-0.0354-1.0E-9 and xp[1]<-0.0202)
        y = xp[1]+0.0354;
      else
      {
        std::cout << xp[1] << std::endl;
        //dserror("coordinates do not match to ORACLES problem");
      }
      // compute y+
      double yplus = y/ltau;

      // generate noise via random number between -1 and +1
      randomnumber = DRT::Problem::Instance()->Random()->Uni();
      noise = fluc * randomnumber;

      // return velocity value in x-direction
      if (yplus <= 5.0)
        return utau * yplus + noise;
      else if ((yplus > 5.0) and (yplus <= 30.0))
        return utau * ( 5.0 * log(yplus) - 3.05 ) + noise;
      else if (yplus > 30.0)
      {
        double upre = utau * ( 2.5 * log(yplus) + 5.5 );
        if (upre > umax) return umax + noise;
        else             return upre + noise;
      }
    }
    // return velocity value in z or y-direction
    case 1:
    case 2:
    {
      // generate noise via random number
      randomnumber = DRT::Problem::Instance()->Random()->Uni();
      noise = umax * fluc * randomnumber;

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
double FLD::JefferyHamelFlowFunction::RadialVelocity(
    const double& theta
    )
{
#ifdef DEBUG
  if (theta < -0.01 or theta > (M_PI/4.0+0.01))
  {
    dserror("angle out of range! Must be between 0 and PI/4");
  }
#endif

  const double alpha = theta - (M_PI/8.0);

  // generated by Mathematica 6
  return -170.00000000402173
        + 51.73882280936082*pow(alpha, 2)
        + 1448.6788580912125*pow(alpha, 4)
        + 16137.997882926256*pow(alpha, 6)
        + 93895.51525051043*pow(alpha, 8)
        + 327716.92421506916*pow(alpha,10)
        - 507543.9679454239*pow(alpha,12)
        + 2.5804347181453653e7*pow(alpha,14)
        - 6.684741949638646e8* pow(alpha,16)
        + 1.0642871526439642e10*pow(alpha,18)
        - 1.2944120060793152e11*pow(alpha,20)
        + 1.1321070844277632e12*pow(alpha,22)
        - 7.135693178454414e12* pow(alpha,24)
        + 3.151583189780962e13*pow(alpha,26)
        - 9.234803801165095e13* pow(alpha,28)
        + 1.6246696343829447e14*pow(alpha,30)
        - 1.3098985134542753e14* pow(alpha,32);
}


double FLD::JefferyHamelFlowFunction::Evaluate(const int index, const double* xp,
    double)
{

  const double x = xp[0];
  const double y = xp[1];

  const double theta = atan(y/x);
  const double u_theta = RadialVelocity(theta);

  // you need to multiply with your viscosity outside of this function to
  // generate flow with higher Reynolds-numbers
  const double nu = 1;

  switch (index)
  {
  case 0:
    return  nu * (u_theta/(x*x+y*y))*x;
  case 1:
    return  nu * (u_theta/(x*x+y*y))*y;
  case 2:
    return 0.0;
  case 3:
    return 0.0;
  default:
    return 0.0;
  }
}
