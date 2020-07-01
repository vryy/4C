/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of functions for structure problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "str_functions.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/stvenantkirchhoff.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STR::WeaklyCompressibleEtienneFSIStructureFunction::WeaklyCompressibleEtienneFSIStructureFunction(
    int mat_id_struc)
    : Function(), youngmodulus_(0.0), poissonratio_(0.0), strucdensity_(0.0)
{
  // get materials
  Teuchos::RCP<MAT::PAR::Material> mat_struc =
      DRT::Problem::Instance()->Materials()->ById(mat_id_struc);
  if (mat_struc->Type() != INPAR::MAT::m_stvenant)
    dserror("Material %d is not a St.Venant-Kirchhoff structure", mat_id_struc);
  MAT::PAR::Parameter* params_struc = mat_struc->Parameter();
  MAT::PAR::StVenantKirchhoff* fparams_struc =
      dynamic_cast<MAT::PAR::StVenantKirchhoff*>(params_struc);
  if (!fparams_struc) dserror("Material does not cast to St.Venant-Kirchhoff structure");


  // get data
  youngmodulus_ = fparams_struc->youngs_;
  poissonratio_ = fparams_struc->poissonratio_;
  strucdensity_ = fparams_struc->density_;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double STR::WeaklyCompressibleEtienneFSIStructureFunction::Evaluate(
    int index, const double* xp, double t)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];

  // initialize variables
  LINALG::Matrix<2, 1> u_ex;

  // evaluate variables
  u_ex(0) = -((cos(2. * PI * t) * cos(2. * PI * x)) / 6. + 1.) * (y - 1.);
  u_ex(1) = -(sin(2. * PI * x) * sin(2. * PI * (t + 1. / 4.)) * (cos(2. * PI * x) - 1.)) / 20.;

  switch (index)
  {
    case 0:
      return u_ex(0);
    case 1:
      return u_ex(1);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> STR::WeaklyCompressibleEtienneFSIStructureFunction::EvaluateTimeDerivative(
    const int index, const double* xp, const double t, const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(index, xp, t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // ease notation
    double x = xp[0];
    double y = xp[1];

    // initialize variables
    LINALG::Matrix<2, 1> dudt_ex;

    // evaluate variables
    dudt_ex(0) = (PI * cos(2. * PI * x) * sin(2. * PI * t) * (y - 1.)) / 3.;
    dudt_ex(1) =
        -(PI * sin(2. * PI * x) * cos(2. * PI * (t + 1. / 4.)) * (cos(2. * PI * x) - 1.)) / 10.;

    switch (index)
    {
      case 0:
        res[1] = dudt_ex(0);
        break;
      case 1:
        res[1] = dudt_ex(1);
        break;

      default:
        res[1] = 0.0;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
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
STR::WeaklyCompressibleEtienneFSIStructureForceFunction::
    WeaklyCompressibleEtienneFSIStructureForceFunction(int mat_id_struc)
    : Function(), youngmodulus_(0.0), poissonratio_(0.0), strucdensity_(0.0)
{
  // get materials
  Teuchos::RCP<MAT::PAR::Material> mat_struc =
      DRT::Problem::Instance()->Materials()->ById(mat_id_struc);
  if (mat_struc->Type() != INPAR::MAT::m_stvenant)
    dserror("Material %d is not a St.Venant-Kirchhoff structure", mat_id_struc);
  MAT::PAR::Parameter* params_struc = mat_struc->Parameter();
  MAT::PAR::StVenantKirchhoff* fparams_struc =
      dynamic_cast<MAT::PAR::StVenantKirchhoff*>(params_struc);
  if (!fparams_struc) dserror("Material does not cast to St.Venant-Kirchhoff structure");

  // get data
  youngmodulus_ = fparams_struc->youngs_;
  poissonratio_ = fparams_struc->poissonratio_;
  strucdensity_ = fparams_struc->density_;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double STR::WeaklyCompressibleEtienneFSIStructureForceFunction::Evaluate(
    int index, const double* xp, double t)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double E = youngmodulus_;
  double v = poissonratio_;
  double r = strucdensity_;

  // initialize variables
  LINALG::Matrix<2, 1> f_u_ex;

  // evaluate variables
  f_u_ex(0) = (2. * (std::pow(PI, 2.)) * cos(2. * PI * t) * cos(2. * PI * x) * (y - 1.) *
                  (E - r - E * v + r * v + 2. * r * (std::pow(v, 2.)))) /
              (3. * (2. * (std::pow(v, 2.)) + v - 1.));
  f_u_ex(1) = ((std::pow(PI, 2.)) * r * sin(2. * PI * x) * sin(2. * PI * (t + 1. / 4.)) *
                  (cos(2. * PI * x) - 1.)) /
                  5. -
              (E * ((PI * sin(2. * PI * x) * sin(2. * PI * (t + 1. / 4.))) / 3. +
                       (3. * (std::pow(PI, 2.)) * cos(2. * PI * x) * sin(2. * PI * x) *
                           sin(2. * PI * (t + 1. / 4.))) /
                           5. +
                       ((std::pow(PI, 2.)) * sin(2. * PI * x) * sin(2. * PI * (t + 1. / 4.)) *
                           (cos(2. * PI * x) - 1.)) /
                           5.)) /
                  (2. * v + 2.) +
              (E * PI * v * sin(2. * PI * x) * sin(2. * PI * (t + 1. / 4.))) /
                  (3. * (2. * v - 1.) * (v + 1.));

  switch (index)
  {
    case 0:
      return f_u_ex(0);
    case 1:
      return f_u_ex(1);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> STR::WeaklyCompressibleEtienneFSIStructureForceFunction::EvaluateTimeDerivative(
    const int index, const double* xp, const double t, const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(index, xp, t);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}
