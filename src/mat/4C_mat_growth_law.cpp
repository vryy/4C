/*----------------------------------------------------------------------*/
/*! \file
\brief growth law for the volumetric growth materials. This file contains routines needed for the
calculation of the volumetric growth parameter theta. These so called growthlaws are to be used for
the isovolumetric split growth materials Mat::Growth. Upon request a growthlaw delivers a growth
factor \f$\vartheta\f$ and its derivative wrt. \f$\frac{\partial \vartheta}{\partial C}\f$.

\level 2

*/

/*----------------------------------------------------------------------*/


#include "4C_mat_growth_law.hpp"

#include "4C_comm_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_growth.hpp"
#include "4C_mat_service.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
Mat::GrowthLaw::GrowthLaw() : params_(nullptr) {}

/*----------------------------------------------------------------------------*/
Mat::GrowthLaw::GrowthLaw(Core::Mat::PAR::Parameter* params) : params_(params) {}

/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawDyn::GrowthLawDyn(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata), abstol_(matdata.parameters.Get<double>("TOL"))
{
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::GrowthLawDyn::create_material()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawDyn::GrowthLawDyn() : GrowthLaw() {}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawDyn::GrowthLawDyn(Core::Mat::PAR::Parameter* params) : GrowthLaw(params) {}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawDyn::Evaluate(double* thetainit, const double& thetaold,
    Core::LinAlg::Matrix<6, 1>* dthetadC, Mat::Growth& matgrowth,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<6, 1>* glstrain,
    const Core::LinAlg::Matrix<3, 1>& refdir, const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // get gauss point number
  double dt = params.get<double>("delta time", -1.0);

  double theta = *thetainit;

  Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatdach(true);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Sdachvec(true);

  // identity tensor I
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Cvec(*glstrain);
  Cvec.Scale(2.0);
  Cvec += Id;

  // calculate growth part F_g of the deformation gradient F
  Core::LinAlg::Matrix<3, 3> F_g(true);
  CalcFg(theta, thetaold, gp, defgrd, refdir, curdir, histdefgrd, F_g);
  // calculate F_g^(-1)
  Core::LinAlg::Matrix<3, 3> F_ginv(true);
  F_ginv.Invert(F_g);
  // elastic deformation gradient F_e = F * F_g^(-1)
  Core::LinAlg::Matrix<3, 3> defgrddach(true);  //*defgrd);
  defgrddach.MultiplyNN(*defgrd, F_ginv);       // Scale(1.0 / theta);
  // elastic right Cauchy-Green Tensor Cdach = F_e^T * F_e (= F_g^-T C F_g^-1)
  Core::LinAlg::Matrix<3, 3> Cdach(true);
  Cdach.MultiplyTN(defgrddach, defgrddach);

  // transform Cdach into a vector
  Core::LinAlg::Matrix<6, 1> Cdachvec(true);
  Core::LinAlg::Voigt::Strains::matrix_to_vector(Cdach, Cdachvec);

  // elastic Green Lagrange strain
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> glstraindachvec(Cdachvec);
  glstraindachvec -= Id;
  glstraindachvec.Scale(0.5);
  // elastic 2 PK stress and constitutive matrix
  matgrowth.EvaluateElastic(
      &defgrddach, &glstraindachvec, &Sdachvec, &cmatdach, params, gp, eleGID);

  // the growth trigger (i.e. stress, strain, ...)
  double growthtrigger = 0.0;
  evaluate_growth_trigger(growthtrigger, Sdachvec, Cdachvec, refdir, theta, consttrig);

  // evaluate growth function
  double growthfunc = 0.0;
  double dgrowthfunctheta = 0.0;
  evaluate_growth_function(growthfunc, growthtrigger, theta);
  evaluate_growth_function_deriv_theta(
      dgrowthfunctheta, growthtrigger, theta, Cdachvec, Sdachvec, cmatdach, F_g, refdir);

  double residual = thetaold - theta + growthfunc * dt;

  int localistep = 0;
  double thetaquer = 0.0;
  int maxstep = 30;

  // local Newton iteration to obtain exact theta
  while (abs(residual) > Parameter()->abstol_ && localistep < maxstep)
  {
    localistep += 1;

    // evaluate derivative of growth function w.r.t. growth factor
    thetaquer = 1.0 - dgrowthfunctheta * dt;

    // damping strategy
    double omega = 2.0;
    double thetatemp = theta;
    double residualtemp = residual;
    double omegamin = 1.0 / 64.0;
    while (abs(residualtemp) > (1.0 - 0.5 * omega) * abs(residual) and omega > omegamin)
    {
      // update of theta
      omega = 0.5 * omega;
      thetatemp = theta + omega * residual / thetaquer;

      // calculate growth part F_g of the deformation gradient F
      F_g.putScalar(0.0);
      CalcFg(thetatemp, thetaold, gp, defgrd, refdir, curdir, histdefgrd, F_g);
      // calculate F_g^(-1)
      F_ginv.putScalar(0.0);
      F_ginv.Invert(F_g);
      // elastic deformation gradient F_e = F * F_g^(-1)
      defgrddach.putScalar(0.0);
      defgrddach.MultiplyNN(*defgrd, F_ginv);
      // elastic right Cauchy-Green tensor Cdach = F_e^T * F_e (= F_g^-T C F_g^-1)
      Cdach.putScalar(0.0);
      Cdach.MultiplyTN(defgrddach, defgrddach);

      // transform Cdach into a vector
      Cdachvec.putScalar(0.0);
      Core::LinAlg::Voigt::Strains::matrix_to_vector(Cdach, Cdachvec);

      glstraindachvec = Cdachvec;
      glstraindachvec -= Id;
      glstraindachvec.Scale(0.5);
      cmatdach.putScalar(0.0);
      Sdachvec.putScalar(0.0);
      matgrowth.EvaluateElastic(
          &defgrddach, &glstraindachvec, &Sdachvec, &cmatdach, params, gp, eleGID);

      // the growth trigger (i.e. stress, strain, ...)
      growthtrigger = 0.0;
      evaluate_growth_trigger(growthtrigger, Sdachvec, Cdachvec, refdir, thetatemp, consttrig);

      growthfunc = 0.0;
      evaluate_growth_function(growthfunc, growthtrigger, thetatemp);

      residualtemp = thetaold - thetatemp + growthfunc * dt;
    }  // end of damping loop
    residual = residualtemp;
    theta = thetatemp;

    // evaluate derivative of growth function w.r.t. growth factor
    evaluate_growth_function_deriv_theta(
        dgrowthfunctheta, growthtrigger, theta, Cdachvec, Sdachvec, cmatdach, F_g, refdir);

    if (omega <= omegamin and abs(residualtemp) > (1.0 - 0.5 * omega) * abs(residual))
    {
      std::cout << gp << ": Theta " << thetatemp << " residual " << residualtemp << " stress "
                << growthtrigger << std::endl;
      // FOUR_C_THROW("no damping coefficient found");
    }

  }  // end of local Newton iteration


  if (localistep == maxstep && abs(residual) > Parameter()->abstol_)
  {
    FOUR_C_THROW(
        "local Newton iteration did not converge after %i steps: residual: %e, thetaold: %f,"
        " theta:  %f, growthtrigger: %e",
        maxstep, residual, thetaold, theta, growthtrigger);
  }

  thetaquer = 1.0 - dgrowthfunctheta * dt;

  // constitutive matrix including growth cmat = F_g^-1 F_g^-1 cmatdach F_g^-T F_g^-T
  Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelastic(true);
  cmatelastic = Mat::PullBackFourTensor<3>(F_ginv, cmatdach);
  // calculate stress
  // 2PK stress S = F_g^-1 Sdach F_g^-T
  Core::LinAlg::Matrix<3, 3> Sdach(true);
  Core::LinAlg::Voigt::Stresses::vector_to_matrix(Sdachvec, Sdach);

  Core::LinAlg::Matrix<3, 3> tmp(true);
  tmp.MultiplyNT(Sdach, F_ginv);
  Core::LinAlg::Matrix<3, 3> S(true);
  S.MultiplyNN(F_ginv, tmp);

  Core::LinAlg::Matrix<6, 1> Svec(true);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(S, Svec);

  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> dgrowthfuncdCvec(true);
  evaluate_growth_function_deriv_c(
      dgrowthfuncdCvec, growthtrigger, theta, Cvec, Svec, cmatelastic, refdir);

  // d theta/ d C
  dthetadC->Update(dt / thetaquer, dgrowthfuncdCvec);

  *thetainit = theta;
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawDyn::EvaluatePDeriv(double* thetainit, const double& thetaold,
    Teuchos::RCP<Mat::So3Material> matelastic, const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params, const int eleGID)
{
  FOUR_C_THROW("This must still be implemented for dynamic type growth laws.");
}

/*----------------------------------------------------------------------------*/
double Mat::GrowthLawDyn::DensityScale(const double theta)
{
  // nothing to do
  return 0;
}
double Mat::GrowthLawDyn::DensityDerivScale(const double theta)
{
  // nothing to do
  return 0;
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawStatic::GrowthLawStatic() : GrowthLaw() {}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawStatic::GrowthLawStatic(Core::Mat::PAR::Parameter* params) : GrowthLaw(params) {}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawStatic::EvaluatePDeriv(double* thetainit, const double& thetaold,
    Teuchos::RCP<Mat::So3Material> matelastic, const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params, const int eleGID)
{
  FOUR_C_THROW("Use implementation of a specific growth law overriding this function.");
}

/*----------------------------------------------------------------------------*/
double Mat::GrowthLawStatic::DensityScale(const double theta)
{
  // nothing to do
  return 0;
}
/*----------------------------------------------------------------------------*/
double Mat::GrowthLawStatic::DensityDerivScale(const double theta)
{
  // nothing to do
  return 0;
}



/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::GrowthLawAnisoStrain::create_material()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Mat::GrowthLaw> Mat::PAR::GrowthLawAnisoStrain::CreateGrowthLaw()
{
  return Teuchos::rcp(new Mat::GrowthLawAnisoStrain(this));
}

/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawAnisoStrain::GrowthLawAnisoStrain(const Core::Mat::PAR::Parameter::Data& matdata)
    : GrowthLawDyn(matdata),
      tau_(matdata.parameters.Get<double>("TAU")),
      taurev_(matdata.parameters.Get<double>("TAU_REV")),
      thetamin_(matdata.parameters.Get<double>("THETA_MIN")),
      thetamax_(matdata.parameters.Get<double>("THETA_MAX")),
      gamma_(matdata.parameters.Get<double>("GAMMA")),
      gammarev_(matdata.parameters.Get<double>("GAMMA_REV")),
      lambdacrit_(matdata.parameters.Get<double>("LAMBDA_CRIT"))
{
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawAnisoStrain::GrowthLawAnisoStrain() : GrowthLawDyn(nullptr) {}


/*----------------------------------------------------------------------------*/
Mat::GrowthLawAnisoStrain::GrowthLawAnisoStrain(Mat::PAR::GrowthLawAnisoStrain* params)
    : GrowthLawDyn(params)
{
}



/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStrain::evaluate_growth_function(
    double& growthfunc, const double growthtrig, const double theta)
{
  Mat::PAR::GrowthLawAnisoStrain* params = Parameter();
  // parameters
  const double tau = params->tau_;
  const double taurev = params->taurev_;
  const double thetamin = params->thetamin_;
  const double thetamax = params->thetamax_;
  const double gamma = params->gamma_;
  const double gammarev = params->gammarev_;
  const double lambdacrit = params->lambdacrit_;

  double ktheta = 0.0;

  if (growthtrig >= lambdacrit)
    ktheta = (1. / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma);
  else
    ktheta = (1. / taurev) * pow(((theta - thetamin) / (thetamax - thetamin)), gammarev);

  growthfunc = ktheta * (growthtrig - lambdacrit);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStrain::evaluate_growth_trigger(double& growthtrig,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<6, 1>& Cdachvec, const Core::LinAlg::Matrix<3, 1>& direction,
    const double& theta, const double& consttrig)
{
  // transform Cdachvec into a matrix
  Core::LinAlg::Matrix<3, 3> Cdach(true);
  Core::LinAlg::Voigt::Strains::vector_to_matrix(Cdachvec, Cdach);

  Core::LinAlg::Matrix<3, 1> CdachDir(true);
  CdachDir.MultiplyNN(1.0, Cdach, direction);

  // elastic fiber stretch
  growthtrig =
      sqrt(CdachDir(0) * direction(0) + CdachDir(1) * direction(1) + CdachDir(2) * direction(2));
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStrain::evaluate_growth_function_deriv_theta(double& dgrowthfunctheta,
    double growthtrig, double theta, const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const Core::LinAlg::Matrix<3, 3>& F_g, const Core::LinAlg::Matrix<3, 1>& direction)
{
  Mat::PAR::GrowthLawAnisoStrain* params = Parameter();
  // parameters
  const double tau = params->tau_;
  const double taurev = params->taurev_;
  const double thetamin = params->thetamin_;
  const double thetamax = params->thetamax_;
  const double gamma = params->gamma_;
  const double gammarev = params->gammarev_;
  const double lambdacrit = params->lambdacrit_;

  double ktheta = 0.0;
  double dktheta = 0.0;

  if (growthtrig >= lambdacrit)
  {
    ktheta = (1. / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma);
    dktheta = -(gamma / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma - 1.) /
              (thetamax - thetamin);
  }
  else
  {
    ktheta = (1. / taurev) * pow(((theta - thetamin) / (thetamax - thetamin)), gammarev);
    dktheta = (gammarev / taurev) *
              pow(((theta - thetamin) / (thetamax - thetamin)), gammarev - 1.) /
              (thetamax - thetamin);
  }

  dgrowthfunctheta = dktheta * (growthtrig - lambdacrit) + ktheta * (growthtrig / theta);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStrain::evaluate_growth_function_deriv_c(
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double growthtrig, double theta,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Svec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
    const Core::LinAlg::Matrix<3, 1>& direction)
{
  Mat::PAR::GrowthLawAnisoStrain* params = Parameter();
  // parameters
  const double tau = params->tau_;
  const double taurev = params->taurev_;
  const double thetamin = params->thetamin_;
  const double thetamax = params->thetamax_;
  const double gamma = params->gamma_;
  const double gammarev = params->gammarev_;
  const double lambdacrit = params->lambdacrit_;

  double ktheta = 0.0;

  if (growthtrig >= lambdacrit)
    ktheta = (1. / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma);
  else
    ktheta = (1. / taurev) * pow(((theta - thetamin) / (thetamax - thetamin)), gammarev);


  Core::LinAlg::Matrix<3, 3> dgrowthfuncdC(true);
  dgrowthfuncdC.MultiplyNT(direction, direction);

  dgrowthfuncdC.Scale(ktheta / (theta * theta * 2. * growthtrig));

  // transform dgrowthfuncdC into a vector in Voigt notation
  dgrowthfuncdCvec(0) = dgrowthfuncdC(0, 0);
  dgrowthfuncdCvec(1) = dgrowthfuncdC(1, 1);
  dgrowthfuncdCvec(2) = dgrowthfuncdC(2, 2);
  dgrowthfuncdCvec(3) = dgrowthfuncdC(0, 1);
  dgrowthfuncdCvec(4) = dgrowthfuncdC(1, 2);
  dgrowthfuncdCvec(5) = dgrowthfuncdC(0, 2);
}

///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStrain::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<3, 1>& refdir,
    const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, Core::LinAlg::Matrix<3, 3>& F_g)
{
  Core::LinAlg::Matrix<3, 3> fdf(true);
  fdf.MultiplyNT(refdir, refdir);

  for (int i = 0; i < 3; i++) F_g(i, i) = 1.0;
  F_g.Update(theta - 1.0, fdf, 1.0);
}


/*----------------------------------------------------------------------------*/
double Mat::GrowthLawAnisoStrain::DensityScale(const double theta)
{
  // in anisotropic uni-directional growth case: rho_0 = theta * rho_dach,
  // thus d(rho_0)/d(theta) = rho_dach
  return theta;
}
/*----------------------------------------------------------------------------*/
double Mat::GrowthLawAnisoStrain::DensityDerivScale(const double theta) { return 1.0; }


/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::GrowthLawAnisoStress::create_material()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Mat::GrowthLaw> Mat::PAR::GrowthLawAnisoStress::CreateGrowthLaw()
{
  return Teuchos::rcp(new Mat::GrowthLawAnisoStress(this));
}

/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawAnisoStress::GrowthLawAnisoStress(const Core::Mat::PAR::Parameter::Data& matdata)
    : GrowthLawDyn(matdata),
      tau_(matdata.parameters.Get<double>("TAU")),
      taurev_(matdata.parameters.Get<double>("TAU_REV")),
      thetamin_(matdata.parameters.Get<double>("THETA_MIN")),
      thetamax_(matdata.parameters.Get<double>("THETA_MAX")),
      gamma_(matdata.parameters.Get<double>("GAMMA")),
      gammarev_(matdata.parameters.Get<double>("GAMMA_REV")),
      pcrit_(matdata.parameters.Get<double>("P_CRIT"))
{
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawAnisoStress::GrowthLawAnisoStress() : GrowthLawDyn(nullptr) {}


/*----------------------------------------------------------------------------*/
Mat::GrowthLawAnisoStress::GrowthLawAnisoStress(Mat::PAR::GrowthLawAnisoStress* params)
    : GrowthLawDyn(params)
{
}



/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStress::evaluate_growth_function(
    double& growthfunc, const double growthtrig, const double theta)
{
  Mat::PAR::GrowthLawAnisoStress* params = Parameter();
  // parameters
  const double tau = params->tau_;
  const double taurev = params->taurev_;
  const double thetamin = params->thetamin_;
  const double thetamax = params->thetamax_;
  const double gamma = params->gamma_;
  const double gammarev = params->gammarev_;
  const double pcrit = params->pcrit_;

  double ktheta = 0.0;
  if (growthtrig >= pcrit)
    ktheta = (1. / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma);
  else
    ktheta = (1. / taurev) * pow(((theta - thetamin) / (thetamax - thetamin)), gammarev);

  growthfunc = ktheta * (growthtrig - pcrit);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStress::evaluate_growth_trigger(double& growthtrig,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<6, 1>& Cdachvec, const Core::LinAlg::Matrix<3, 1>& direction,
    const double& theta, const double& consttrig)
{
  // trace of elastic Mandel stress tensor
  growthtrig = Cdachvec(0) * Sdachvec(0) + Cdachvec(1) * Sdachvec(1) + Cdachvec(2) * Sdachvec(2) +
               Cdachvec(3) * Sdachvec(3) + Cdachvec(4) * Sdachvec(4) + Cdachvec(5) * Sdachvec(5);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStress::evaluate_growth_function_deriv_theta(double& dgrowthfunctheta,
    double growthtrig, double theta, const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const Core::LinAlg::Matrix<3, 3>& F_g, const Core::LinAlg::Matrix<3, 1>& direction)
{
  Mat::PAR::GrowthLawAnisoStress* params = Parameter();
  // parameters
  const double tau = params->tau_;
  const double taurev = params->taurev_;
  const double thetamin = params->thetamin_;
  const double thetamax = params->thetamax_;
  const double gamma = params->gamma_;
  const double gammarev = params->gammarev_;
  const double pcrit = params->pcrit_;

  double ktheta = 0.0;
  double dktheta = 0.0;

  if (growthtrig > pcrit)
  {
    ktheta = (1. / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma);
    dktheta = -(gamma / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma - 1.) /
              (thetamax - thetamin);
  }
  else
  {
    ktheta = (1. / taurev) * pow(((theta - thetamin) / (thetamax - thetamin)), gammarev);
    dktheta = (gammarev / taurev) *
              pow(((theta - thetamin) / (thetamax - thetamin)), gammarev - 1.) /
              (thetamax - thetamin);
  }

  // calculate F_g^(-1)
  Core::LinAlg::Matrix<3, 3> F_ginv(true);
  F_ginv.Invert(F_g);

  // calculate d(F_g)/d(theta) = 1 - f_0 \otimes f_0
  Core::LinAlg::Matrix<3, 3> fdf(true);
  fdf.MultiplyNT(direction, direction);

  Core::LinAlg::Matrix<3, 3> dFgdtheta(true);
  for (int i = 0; i < 3; i++) dFgdtheta(i, i) = 1.0;

  dFgdtheta.Update(-1.0, fdf, 1.0);


  // transform Cdachvec into a matrix
  Core::LinAlg::Matrix<3, 3> Cdach(true);
  Core::LinAlg::Voigt::Strains::vector_to_matrix(Cdachvec, Cdach);

  Core::LinAlg::Matrix<3, 3> dFgT_Cdach(true);
  dFgT_Cdach.MultiplyTN(dFgdtheta, Cdach);
  Core::LinAlg::Matrix<3, 3> dCdachdtheta(true);
  dCdachdtheta.MultiplyTN(F_ginv, dFgT_Cdach);

  Core::LinAlg::Matrix<3, 3> dFg_Fginv(true);
  dFg_Fginv.MultiplyNN(dFgdtheta, F_ginv);
  Core::LinAlg::Matrix<3, 3> Cdach_dFg_Fginv(true);
  Cdach_dFg_Fginv.MultiplyNN(Cdach, dFg_Fginv);

  dCdachdtheta.Update(-1.0, Cdach_dFg_Fginv, -1.0);

  // transform dCdachdtheta into a vector
  Core::LinAlg::Matrix<6, 1> dCdachdthetavec(true);
  Core::LinAlg::Voigt::Strains::matrix_to_vector(dCdachdtheta, dCdachdthetavec);

  Core::LinAlg::Matrix<6, 1> dSdachdthetavec(true);

  for (int i = 0; i < 6; i++)
  {
    dSdachdthetavec(i) =
        0.5 * (cmatelastic(i, 0) * dCdachdthetavec(0) + cmatelastic(i, 1) * dCdachdthetavec(1) +
                  cmatelastic(i, 2) * dCdachdthetavec(2) + cmatelastic(i, 3) * dCdachdthetavec(3) +
                  cmatelastic(i, 4) * dCdachdthetavec(4) + cmatelastic(i, 5) * dCdachdthetavec(5));
  }

  // derivative of trace of Mandel stress tensor w.r.t. theta
  // cf. Goektepe et al. J Theor Biol (2010), eq. 26 ff.
  const double dgrowthtrig = dCdachdthetavec(0) * Sdachvec(0) + dCdachdthetavec(1) * Sdachvec(1) +
                             dCdachdthetavec(2) * Sdachvec(2) + dCdachdthetavec(3) * Sdachvec(3) +
                             dCdachdthetavec(4) * Sdachvec(4) + dCdachdthetavec(5) * Sdachvec(5)

                             + dSdachdthetavec(0) * Cdachvec(0) + dSdachdthetavec(1) * Cdachvec(1) +
                             dSdachdthetavec(2) * Cdachvec(2) + dSdachdthetavec(3) * Cdachvec(3) +
                             dSdachdthetavec(4) * Cdachvec(4) + dSdachdthetavec(5) * Cdachvec(5);


  dgrowthfunctheta = dktheta * (growthtrig - pcrit) + ktheta * dgrowthtrig;
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStress::evaluate_growth_function_deriv_c(
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double growthtrig, double theta,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Svec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
    const Core::LinAlg::Matrix<3, 1>& direction)
{
  Mat::PAR::GrowthLawAnisoStress* params = Parameter();
  // parameters
  const double tau = params->tau_;
  const double taurev = params->taurev_;
  const double thetamin = params->thetamin_;
  const double thetamax = params->thetamax_;
  const double gamma = params->gamma_;
  const double gammarev = params->gammarev_;
  const double pcrit = params->pcrit_;

  double ktheta = 0.;

  if (growthtrig >= pcrit)
    ktheta = (1. / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma);
  else
    ktheta = (1. / taurev) * pow(((theta - thetamin) / (thetamax - thetamin)), gammarev);


  for (int j = 0; j < NUM_STRESS_3D; j++)
  {
    double Ccmatelasj = cmat(0, j) * Cvec(0) + cmat(1, j) * Cvec(1) + cmat(2, j) * Cvec(2) +
                        cmat(3, j) * Cvec(3) + cmat(4, j) * Cvec(4) + cmat(5, j) * Cvec(5);

    dgrowthfuncdCvec(j) = ktheta * (Svec(j) + 0.5 * Ccmatelasj);
  }
}


///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStress::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<3, 1>& refdir,
    const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, Core::LinAlg::Matrix<3, 3>& F_g)
{
  Core::LinAlg::Matrix<3, 3> fdf(true);
  fdf.MultiplyNT(refdir, refdir);

  for (int i = 0; i < 3; i++) F_g(i, i) = theta;
  F_g.Update(1.0 - theta, fdf, 1.0);
}


/*----------------------------------------------------------------------------*/
double Mat::GrowthLawAnisoStress::DensityScale(const double theta)
{
  // in anisotropic bi-directional growth case: rho_0 = theta^2 * rho_dach,
  // thus d(rho_0)/d(theta) = 2 * theta * rho_dach
  return theta * theta;
}
/*----------------------------------------------------------------------------*/
double Mat::GrowthLawAnisoStress::DensityDerivScale(const double theta) { return 2.0 * theta; }



/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::GrowthLawAnisoStrainConstTrig::create_material()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Mat::GrowthLaw> Mat::PAR::GrowthLawAnisoStrainConstTrig::CreateGrowthLaw()
{
  return Teuchos::rcp(new Mat::GrowthLawAnisoStrainConstTrig(this));
}

/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawAnisoStrainConstTrig::GrowthLawAnisoStrainConstTrig(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : GrowthLawAnisoStrain(matdata)
{
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawAnisoStrainConstTrig::GrowthLawAnisoStrainConstTrig() : GrowthLawAnisoStrain(nullptr)
{
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawAnisoStrainConstTrig::GrowthLawAnisoStrainConstTrig(
    Mat::PAR::GrowthLawAnisoStrainConstTrig* params)
    : GrowthLawAnisoStrain(params)
{
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStrainConstTrig::evaluate_growth_trigger(double& growthtrig,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<6, 1>& Cdachvec, const Core::LinAlg::Matrix<3, 1>& direction,
    const double& theta, const double& consttrig)
{
  growthtrig = consttrig;
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStrainConstTrig::evaluate_growth_function_deriv_theta(
    double& dgrowthfunctheta, double growthtrig, double theta,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const Core::LinAlg::Matrix<3, 3>& F_g, const Core::LinAlg::Matrix<3, 1>& direction)
{
  Mat::PAR::GrowthLawAnisoStrainConstTrig* params = Parameter();
  // parameters
  const double tau = params->tau_;
  const double taurev = params->taurev_;
  const double thetamin = params->thetamin_;
  const double thetamax = params->thetamax_;
  const double gamma = params->gamma_;
  const double gammarev = params->gammarev_;
  const double lambdacrit = params->lambdacrit_;

  double dktheta = 0.0;

  if (growthtrig >= lambdacrit)
  {
    dktheta = -(gamma / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma - 1.) /
              (thetamax - thetamin);
  }
  else
  {
    dktheta = (gammarev / taurev) *
              pow(((theta - thetamin) / (thetamax - thetamin)), gammarev - 1.) /
              (thetamax - thetamin);
  }

  dgrowthfunctheta = dktheta * (growthtrig - lambdacrit);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStrainConstTrig::evaluate_growth_function_deriv_c(
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double growthtrig, double theta,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Svec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
    const Core::LinAlg::Matrix<3, 1>& direction)
{
  for (int j = 0; j < NUM_STRESS_3D; j++)
  {
    dgrowthfuncdCvec(j) = 0.0;
  }
}



/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::GrowthLawAnisoStressConstTrig::create_material()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Mat::GrowthLaw> Mat::PAR::GrowthLawAnisoStressConstTrig::CreateGrowthLaw()
{
  return Teuchos::rcp(new Mat::GrowthLawAnisoStressConstTrig(this));
}

/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawAnisoStressConstTrig::GrowthLawAnisoStressConstTrig(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : GrowthLawAnisoStress(matdata)
{
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawAnisoStressConstTrig::GrowthLawAnisoStressConstTrig() : GrowthLawAnisoStress(nullptr)
{
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawAnisoStressConstTrig::GrowthLawAnisoStressConstTrig(
    Mat::PAR::GrowthLawAnisoStressConstTrig* params)
    : GrowthLawAnisoStress(params)
{
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStressConstTrig::evaluate_growth_trigger(double& growthtrig,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<6, 1>& Cdachvec, const Core::LinAlg::Matrix<3, 1>& direction,
    const double& theta, const double& consttrig)
{
  growthtrig = consttrig;
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStressConstTrig::evaluate_growth_function_deriv_theta(
    double& dgrowthfunctheta, double growthtrig, double theta,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const Core::LinAlg::Matrix<3, 3>& F_g, const Core::LinAlg::Matrix<3, 1>& direction)
{
  Mat::PAR::GrowthLawAnisoStressConstTrig* params = Parameter();
  // parameters
  const double tau = params->tau_;
  const double taurev = params->taurev_;
  const double thetamin = params->thetamin_;
  const double thetamax = params->thetamax_;
  const double gamma = params->gamma_;
  const double gammarev = params->gammarev_;
  const double pcrit = params->pcrit_;

  double dktheta = 0.0;

  if (growthtrig > pcrit)
  {
    dktheta = -(gamma / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma - 1.) /
              (thetamax - thetamin);
  }
  else
  {
    dktheta = (gammarev / taurev) *
              pow(((theta - thetamin) / (thetamax - thetamin)), gammarev - 1.) /
              (thetamax - thetamin);
  }

  dgrowthfunctheta = dktheta * (growthtrig - pcrit);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAnisoStressConstTrig::evaluate_growth_function_deriv_c(
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double growthtrig, double theta,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Svec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
    const Core::LinAlg::Matrix<3, 1>& direction)
{
  for (int j = 0; j < NUM_STRESS_3D; j++)
  {
    dgrowthfuncdCvec(j) = 0.0;
  }
}



/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::GrowthLawIsoStress::create_material()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Mat::GrowthLaw> Mat::PAR::GrowthLawIsoStress::CreateGrowthLaw()
{
  return Teuchos::rcp(new Mat::GrowthLawIsoStress(this));
}

/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawIsoStress::GrowthLawIsoStress(const Core::Mat::PAR::Parameter::Data& matdata)
    : GrowthLawDyn(matdata),
      thetaplus_(matdata.parameters.Get<double>("THETAPLUS")),
      kthetaplus_(matdata.parameters.Get<double>("KPLUS")),
      mthetaplus_(matdata.parameters.Get<double>("MPLUS")),
      thetaminus_(matdata.parameters.Get<double>("THETAMINUS")),
      kthetaminus_(matdata.parameters.Get<double>("KMINUS")),
      mthetaminus_(matdata.parameters.Get<double>("MMINUS")),
      hommandel_(matdata.parameters.Get<double>("HOMMANDEL"))
{
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawIsoStress::GrowthLawIsoStress() : GrowthLawDyn(nullptr) {}


/*----------------------------------------------------------------------------*/
Mat::GrowthLawIsoStress::GrowthLawIsoStress(Mat::PAR::GrowthLawIsoStress* params)
    : GrowthLawDyn(params)
{
}



/*----------------------------------------------------------------------------*/
void Mat::GrowthLawIsoStress::evaluate_growth_function(
    double& growthfunc, const double growthtrig, const double theta)
{
  Mat::PAR::GrowthLawIsoStress* params = Parameter();
  // parameters
  const double hommandel = params->hommandel_;
  const double thetaplus = params->thetaplus_;      // 1.5;
  const double kthetaplus = params->kthetaplus_;    // 1.0; 0.5;
  const double mthetaplus = params->mthetaplus_;    // 2.0; 4.0;
  const double thetaminus = params->thetaminus_;    // 0.5;
  const double kthetaminus = params->kthetaminus_;  // 2.0; 0.25;
  const double mthetaminus = params->mthetaminus_;  // 3.0; 5.0;

  double ktheta = 0.0;

  if (growthtrig >= hommandel)
    ktheta = kthetaplus * pow((thetaplus - theta) / (thetaplus - 1.0), mthetaplus);
  else
    ktheta = kthetaminus * pow((theta - thetaminus) / (1.0 - thetaminus), mthetaminus);


  growthfunc = ktheta * (growthtrig - hommandel);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawIsoStress::evaluate_growth_trigger(double& growthtrig,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<6, 1>& Cdachvec, const Core::LinAlg::Matrix<3, 1>& direction,
    const double& theta, const double& consttrig)
{
  // trace of elastic Mandel stress tensor
  growthtrig = Cdachvec(0) * Sdachvec(0) + Cdachvec(1) * Sdachvec(1) + Cdachvec(2) * Sdachvec(2) +
               Cdachvec(3) * Sdachvec(3) + Cdachvec(4) * Sdachvec(4) + Cdachvec(5) * Sdachvec(5);
}


/*----------------------------------------------------------------------------*/
void Mat::GrowthLawIsoStress::evaluate_growth_function_deriv_theta(double& dgrowthfunctheta,
    double growthtrig, double theta, const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const Core::LinAlg::Matrix<3, 3>& F_g, const Core::LinAlg::Matrix<3, 1>& direction)
{
  Mat::PAR::GrowthLawIsoStress* params = Parameter();
  // parameters
  const double hommandel = params->hommandel_;
  const double thetaplus = params->thetaplus_;
  const double kthetaplus = params->kthetaplus_;
  const double mthetaplus = params->mthetaplus_;
  const double thetaminus = params->thetaminus_;
  const double kthetaminus = params->kthetaminus_;
  const double mthetaminus = params->mthetaminus_;

  double ktheta = 0.0;
  double dktheta = 0.0;

  if (growthtrig >= hommandel)
  {
    ktheta = kthetaplus * pow((thetaplus - theta) / (thetaplus - 1.0), mthetaplus);
    dktheta = mthetaplus * kthetaplus *
              pow((thetaplus - theta) / (thetaplus - 1.0), mthetaplus - 1.0) / (1.0 - thetaplus);
  }
  else
  {
    ktheta = kthetaminus * pow((theta - thetaminus) / (1.0 - thetaminus), mthetaminus);
    dktheta = mthetaminus * kthetaminus *
              pow((theta - thetaminus) / (1.0 - thetaminus), mthetaminus - 1.0) /
              (1.0 - thetaminus);
  }

  double temp = 0.0;
  for (int i = 0; i < 6; i++)
  {
    temp += Cdachvec(i) * (cmatelastic(i, 0) * Cdachvec(0) + cmatelastic(i, 1) * Cdachvec(1) +
                              cmatelastic(i, 2) * Cdachvec(2) + cmatelastic(i, 3) * Cdachvec(3) +
                              cmatelastic(i, 4) * Cdachvec(4) + cmatelastic(i, 5) * Cdachvec(5));
  }

  const double dtraceM = -(2.0 * growthtrig + temp) / theta;

  dgrowthfunctheta = (dktheta * (growthtrig - hommandel) + ktheta * dtraceM);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawIsoStress::evaluate_growth_function_deriv_c(
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double growthtrig, double theta,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Cvec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& Svec,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
    const Core::LinAlg::Matrix<3, 1>& direction)
{
  Mat::PAR::GrowthLawIsoStress* params = Parameter();
  // parameters
  const double hommandel = params->hommandel_;
  const double thetaplus = params->thetaplus_;      // 1.5;
  const double kthetaplus = params->kthetaplus_;    // 1.0; 0.5;
  const double mthetaplus = params->mthetaplus_;    // 2.0; 4.0;
  const double thetaminus = params->thetaminus_;    // 0.5;
  const double kthetaminus = params->kthetaminus_;  // 2.0; 0.25;
  const double mthetaminus = params->mthetaminus_;  // 3.0; 5.0;

  double ktheta = 0.0;

  if (growthtrig > hommandel)
  {
    ktheta = kthetaplus * pow((thetaplus - theta) / (thetaplus - 1.0), mthetaplus);
  }
  else if (growthtrig < hommandel)
  {
    ktheta = kthetaminus * pow((theta - thetaminus) / (1.0 - thetaminus), mthetaminus);
  }

  for (int j = 0; j < NUM_STRESS_3D; j++)
  {
    double Ccmatelasj = cmat(0, j) * Cvec(0) + cmat(1, j) * Cvec(1) + cmat(2, j) * Cvec(2) +
                        cmat(3, j) * Cvec(3) + cmat(4, j) * Cvec(4) + cmat(5, j) * Cvec(5);

    dgrowthfuncdCvec(j) = ktheta * (Svec(j) + 0.5 * Ccmatelasj);
  }
}


///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void Mat::GrowthLawIsoStress::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<3, 1>& refdir,
    const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, Core::LinAlg::Matrix<3, 3>& F_g)
{
  for (int i = 0; i < 3; i++) F_g(i, i) = theta;
}


/*----------------------------------------------------------------------------*/
double Mat::GrowthLawIsoStress::DensityScale(const double theta)
{
  // in isotropic growth case: rho_0 = theta^3 * rho_dach,
  // thus d(rho_0)/d(theta) = 3 * theta^2 * rho_dach
  return theta * theta * theta;
}
/*----------------------------------------------------------------------------*/
double Mat::GrowthLawIsoStress::DensityDerivScale(const double theta)
{
  return 3.0 * theta * theta;
}



/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawAC::GrowthLawAC(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      Sc1_(matdata.parameters.Get<int>("SCALAR1")),
      alpha_(matdata.parameters.Get<double>("ALPHA")),
      Sc2_(matdata.parameters.Get<int>("SCALAR2")),
      beta_(matdata.parameters.Get<double>("BETA"))
{
  if (Sc1_ < 1) FOUR_C_THROW("At least on scalar field must induce growth");
  if (alpha_ < 0) FOUR_C_THROW("The influence of scalar field SCALAR1 to growth can't be negativ");
  if (Sc2_ < 1)
    FOUR_C_THROW("If you choose a second scalar field to induce growth, choose a existing one!!");
  if (beta_ < 0) FOUR_C_THROW("The influence of scalar field SCALAR2 to growth can't be negativ");
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::GrowthLawAC::create_material() { return Teuchos::null; }

Teuchos::RCP<Mat::GrowthLaw> Mat::PAR::GrowthLawAC::CreateGrowthLaw()
{
  return Teuchos::rcp(new Mat::GrowthLawAC(this));
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawAC::GrowthLawAC() : GrowthLawStatic(nullptr) {}


/*----------------------------------------------------------------------------*/
Mat::GrowthLawAC::GrowthLawAC(Mat::PAR::GrowthLawAC* params) : GrowthLawStatic(params) {}



/*----------------------------------------------------------------------------*/
void Mat::GrowthLawAC::Evaluate(double* theta, const double& thetaold,
    Core::LinAlg::Matrix<6, 1>* dthetadC, Mat::Growth& matgrowth,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<6, 1>* glstrain,
    const Core::LinAlg::Matrix<3, 1>& refdir, const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // get pointer to vector containing the scalar values at the Gauss points
  Teuchos::RCP<std::vector<std::vector<double>>> concentrations =
      params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
          Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));

  const int Sc1 = Parameter()->Sc1_;
  const double alpha = Parameter()->alpha_;
  // NOTE: if no second scalar is chosen to induce growth, this points out the first scalar field
  const int Sc2 = Parameter()->Sc2_;
  // NOTE: if no second scalar is chosen to induce growth, beta is zero and hence has no influence
  // on growth
  const double beta = Parameter()->beta_;

  const double deltagrowth =
      (alpha * concentrations->at(gp).at(Sc1 - 1) + beta * concentrations->at(gp).at(Sc2 - 1)) *
      defgrd->Determinant();

  *theta = pow(1.0 + deltagrowth, 0.33333333333333333);


  const double tmp = deltagrowth / 6.0 * pow(1.0 + deltagrowth, -0.66666666666666666);

  Core::LinAlg::Matrix<3, 3> C(true);
  C.MultiplyTN(*defgrd, *defgrd);
  Core::LinAlg::Matrix<3, 3> Cinv(true);
  Cinv.Invert(C);

  // linearization of growth law
  Core::LinAlg::Matrix<3, 3> dThetadC(Cinv);
  dThetadC.Scale(tmp);

  // transform dThetadC into a vector in Voigt notation
  (*dthetadC)(0) = dThetadC(0, 0);
  (*dthetadC)(1) = dThetadC(1, 1);
  (*dthetadC)(2) = dThetadC(2, 2);
  (*dthetadC)(3) = dThetadC(0, 1);
  (*dthetadC)(4) = dThetadC(1, 2);
  (*dthetadC)(5) = dThetadC(0, 2);

  // set ratio for potential linear interpolation between two elastic materials
  double lambda = 1.0 / (1.0 + deltagrowth);
  // linearization of ratio for potential linear interpolation of between two elastic materials
  Cinv.Scale(-0.5 * deltagrowth * pow(1.0 + deltagrowth, -2.0));

  // transform into a vector in Voigt notation
  Teuchos::RCP<Core::LinAlg::Matrix<6, 1>> dlambda_dC =
      Teuchos::rcp(new Core::LinAlg::Matrix<6, 1>(true));
  (*dlambda_dC)(0) = Cinv(0, 0);
  (*dlambda_dC)(1) = Cinv(1, 1);
  (*dlambda_dC)(2) = Cinv(2, 2);
  (*dlambda_dC)(3) = Cinv(0, 1);
  (*dlambda_dC)(4) = Cinv(1, 2);
  (*dlambda_dC)(5) = Cinv(0, 2);

  // save values in parameter list
  params.set<double>("lambda", lambda);
  params.set<Teuchos::RCP<Core::LinAlg::Matrix<6, 1>>>("dlambda_dC", dlambda_dC);
}


///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void Mat::GrowthLawAC::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<3, 1>& refdir,
    const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, Core::LinAlg::Matrix<3, 3>& F_g)
{
  for (int i = 0; i < 3; i++) F_g(i, i) = theta;
}


/*----------------------------------------------------------------------------*/
double Mat::GrowthLawAC::DensityScale(const double theta)
{
  // in isotropic growth case: rho_0 = theta^3 * rho_dach,
  // thus d(rho_0)/d(theta) = 3 * theta^2 * rho_dach
  return theta * theta * theta;
}
/*----------------------------------------------------------------------------*/
double Mat::GrowthLawAC::DensityDerivScale(const double theta) { return 3.0 * theta * theta; }



/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawACRadial::GrowthLawACRadial(const Core::Mat::PAR::Parameter::Data& matdata)
    : GrowthLawAC(matdata)
{
}

Teuchos::RCP<Mat::GrowthLaw> Mat::PAR::GrowthLawACRadial::CreateGrowthLaw()
{
  return Teuchos::rcp(new Mat::GrowthLawACRadial(this));
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawACRadial::GrowthLawACRadial() : GrowthLawStatic(nullptr) {}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawACRadial::GrowthLawACRadial(Mat::PAR::GrowthLawAC* params) : GrowthLawStatic(params)
{
}



/*----------------------------------------------------------------------------*/
void Mat::GrowthLawACRadial::Evaluate(double* theta, const double& thetaold,
    Core::LinAlg::Matrix<6, 1>* dthetadC, Mat::Growth& matgrowth,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<6, 1>* glstrain,
    const Core::LinAlg::Matrix<3, 1>& refdir, const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // get pointer to vector containing the scalar values at the Gauss points
  Teuchos::RCP<std::vector<std::vector<double>>> concentrations =
      params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
          Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));

  const int Sc1 = Parameter()->Sc1_;
  const double alpha = Parameter()->alpha_;
  // NOTE: if no second scalar is chosen to induce growth, this points out the first scalar field
  const int Sc2 = Parameter()->Sc2_;
  // NOTE: if no second scalar is chosen to induce growth, beta is zero and hence has no influence
  // on growth
  const double beta = Parameter()->beta_;

  const double deltagrowth =
      (alpha * concentrations->at(gp).at(Sc1 - 1) + beta * concentrations->at(gp).at(Sc2 - 1)) *
      defgrd->Determinant();

  *theta = 1.0 + deltagrowth;


  const double tmp = 0.5 * deltagrowth;

  Core::LinAlg::Matrix<3, 3> C(true);
  C.MultiplyTN(*defgrd, *defgrd);
  Core::LinAlg::Matrix<3, 3> Cinv(true);
  Cinv.Invert(C);

  // linearization of growth law
  Core::LinAlg::Matrix<3, 3> dThetadC(Cinv);
  dThetadC.Scale(tmp);

  // transform dThetadC into a vector in Voigt notation
  (*dthetadC)(0) = dThetadC(0, 0);
  (*dthetadC)(1) = dThetadC(1, 1);
  (*dthetadC)(2) = dThetadC(2, 2);
  (*dthetadC)(3) = dThetadC(0, 1);
  (*dthetadC)(4) = dThetadC(1, 2);
  (*dthetadC)(5) = dThetadC(0, 2);

  // set ratio for potential linear interpolation between two elastic materials
  const double lambda = 1.0 / *theta;
  // linearization of ratio for potential linear interpolation of between two elastic materials
  Cinv.Scale(-0.5 * deltagrowth * pow(1.0 + deltagrowth, -2.0));

  // transform into a vector in Voigt notation
  Teuchos::RCP<Core::LinAlg::Matrix<6, 1>> dlambda_dC =
      Teuchos::rcp(new Core::LinAlg::Matrix<6, 1>(true));
  (*dlambda_dC)(0) = Cinv(0, 0);
  (*dlambda_dC)(1) = Cinv(1, 1);
  (*dlambda_dC)(2) = Cinv(2, 2);
  (*dlambda_dC)(3) = Cinv(0, 1);
  (*dlambda_dC)(4) = Cinv(1, 2);
  (*dlambda_dC)(5) = Cinv(0, 2);

  // save values in parameter list
  params.set<double>("lambda", lambda);
  params.set<Teuchos::RCP<Core::LinAlg::Matrix<6, 1>>>("dlambda_dC", dlambda_dC);
}



///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void Mat::GrowthLawACRadial::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<3, 1>& refdir,
    const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, Core::LinAlg::Matrix<3, 3>& F_g

)
{
  Core::LinAlg::Matrix<3, 3> ndn(true);
  ndn.MultiplyNT(curdir.at(gp), curdir.at(gp));

  Core::LinAlg::Matrix<3, 3> F_g_incr(true);
  for (int i = 0; i < 3; i++) F_g_incr(i, i) = 1.0;
  F_g_incr.Update((theta - thetaold) / thetaold, ndn, 1.0);

  F_g.MultiplyNN(F_g_incr, histdefgrd.at(gp));
}



/*----------------------------------------------------------------------------*/
double Mat::GrowthLawACRadial::DensityScale(const double theta)
{
  // in anisotropic uni-directional growth case: rho_0 = theta * rho_dach,
  // thus d(rho_0)/d(theta) = rho_dach
  return theta;
}
/*----------------------------------------------------------------------------*/
double Mat::GrowthLawACRadial::DensityDerivScale(const double theta) { return 1.0; }



/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawACRadialRefConc::GrowthLawACRadialRefConc(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : GrowthLawAC(matdata)
{
}

Teuchos::RCP<Mat::GrowthLaw> Mat::PAR::GrowthLawACRadialRefConc::CreateGrowthLaw()
{
  return Teuchos::rcp(new Mat::GrowthLawACRadialRefConc(this));
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawACRadialRefConc::GrowthLawACRadialRefConc() : GrowthLawStatic(nullptr) {}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawACRadialRefConc::GrowthLawACRadialRefConc(Mat::PAR::GrowthLawAC* params)
    : GrowthLawStatic(params)
{
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawACRadialRefConc::Evaluate(double* theta, const double& thetaold,
    Core::LinAlg::Matrix<6, 1>* dthetadC, Mat::Growth& matgrowth,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<6, 1>* glstrain,
    const Core::LinAlg::Matrix<3, 1>& refdir, const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // get pointer to vector containing the scalar values at the Gauss points
  Teuchos::RCP<std::vector<std::vector<double>>> concentrations =
      params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
          Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));

  const int Sc1 = Parameter()->Sc1_;
  const double alpha = Parameter()->alpha_;
  // NOTE: if no second scalar is chosen to induce growth, this points out the first scalar field
  const int Sc2 = Parameter()->Sc2_;
  // NOTE: if no second scalar is chosen to induce growth, beta is zero and hence has no influence
  // on growth
  const double beta = Parameter()->beta_;

  // const double deltagrowth = (alpha*GetCfac().at(Sc1-1) + beta *
  // GetCfac().at(Sc2-1))*defgrd->Determinant();
  const double deltagrowth =
      (alpha * concentrations->at(gp).at(Sc1 - 1) + beta * concentrations->at(gp).at(Sc2 - 1)) *
      1.0;

  *theta = 1.0 + deltagrowth;

  dthetadC->PutScalar(0.0);

  // set ratio for potential linear interpolation between two elastic materials
  const double lambda = 1.0 / *theta;

  // save values in parameter list
  params.set<double>("lambda", lambda);

  // here we don't need this:
  //  params.set< Teuchos::RCP<Core::LinAlg::Matrix<6,1> > >("dlambda_dC",dlambda_dC);
}


///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void Mat::GrowthLawACRadialRefConc::CalcFg(const double& theta, const double& thetaold,
    const int& gp, const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<3, 1>& refdir, const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, Core::LinAlg::Matrix<3, 3>& F_g)
{
  Core::LinAlg::Matrix<3, 3> ndn(true);
  ndn.MultiplyNT(curdir.at(gp), curdir.at(gp));

  Core::LinAlg::Matrix<3, 3> F_g_incr(true);
  for (int i = 0; i < 3; i++) F_g_incr(i, i) = 1.0;
  F_g_incr.Update((theta - thetaold) / thetaold, ndn, 1.0);

  F_g.MultiplyNN(F_g_incr, histdefgrd.at(gp));
}


/*----------------------------------------------------------------------------*/
double Mat::GrowthLawACRadialRefConc::DensityScale(const double theta)
{
  // in anisotropic uni-directional growth case: rho_0 = theta * rho_dach,
  // thus d(rho_0)/d(theta) = rho_dach
  return theta;
}
/*----------------------------------------------------------------------------*/
double Mat::GrowthLawACRadialRefConc::DensityDerivScale(const double theta) { return 1.0; }



/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::GrowthLawConst::create_material()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Mat::GrowthLaw> Mat::PAR::GrowthLawConst::CreateGrowthLaw()
{
  return Teuchos::rcp(new Mat::GrowthLawConst(this));
}

/*----------------------------------------------------------------------------*/
Mat::PAR::GrowthLawConst::GrowthLawConst(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata)
{
  Epetra_Map dummy_map(1, 1, 0, *(Global::Problem::Instance()->GetCommunicators()->LocalComm()));
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(thetarate)->PutScalar(matdata.parameters.Get<double>("THETARATE"));
}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawConst::GrowthLawConst() : GrowthLawStatic(nullptr) {}

/*----------------------------------------------------------------------------*/
Mat::GrowthLawConst::GrowthLawConst(Mat::PAR::GrowthLawConst* params) : GrowthLawStatic(params) {}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawConst::Evaluate(double* theta, const double& thetaold,
    Core::LinAlg::Matrix<6, 1>* dthetadC, Mat::Growth& matgrowth,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<6, 1>* glstrain,
    const Core::LinAlg::Matrix<3, 1>& refdir, const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  double dt = params.get<double>("delta time", -1.0);
  // map in GetParameter can now calculate LID, so we do not need it here       05/2017 birzle
  // int eleID = Global::Problem::Instance()->GetDis("structure")->ElementColMap()->LID(eleGID);
  *theta = thetaold + Parameter()->GetParameter(Parameter()->thetarate, eleGID) * dt;
  //*theta = Parameter()->GetParameter(Parameter()->thetarate,eleID);

  dthetadC->putScalar(0.0);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthLawConst::EvaluatePDeriv(double* theta, const double& thetaold,
    Teuchos::RCP<Mat::So3Material> matelastic, const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params, const int eleGID)
{
  double dt = params.get<double>("delta time", -1.0);

  int deriv = params.get<int>("matparderiv", -1);
  if (deriv == Parameter()->thetarate)
  {
    *theta = dt;
  }
  else if (deriv == -1)
    FOUR_C_THROW("you should only end up here with a valid matparderiv flag!");
}

///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void Mat::GrowthLawConst::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<3, 1>& refdir,
    const std::vector<Core::LinAlg::Matrix<3, 1>>& curdir,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& histdefgrd, Core::LinAlg::Matrix<3, 3>& F_g)
{
  for (int i = 0; i < 3; i++) F_g(i, i) = theta;
}


/*----------------------------------------------------------------------------*/
double Mat::GrowthLawConst::DensityScale(const double theta)
{
  // in isotropic growth case: rho_0 = theta^3 * rho_dach,
  // thus d(rho_0)/d(theta) = 3 * theta^2 * rho_dach
  return theta * theta * theta;
}
/*----------------------------------------------------------------------------*/
double Mat::GrowthLawConst::DensityDerivScale(const double theta) { return 3.0 * theta * theta; }

FOUR_C_NAMESPACE_CLOSE
