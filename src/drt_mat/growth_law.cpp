/*----------------------------------------------------------------------*/
/*! \file
\brief growth law for the volumetric growth materials. This file contains routines needed for the
calculation of the volumetric growth parameter theta. These so called growthlaws are to be used for
the isovolumetric split growth materials MAT::Growth. Upon request a growthlaw delivers a growth
factor \f$\vartheta\f$ and its derivative wrt. \f$\frac{\partial \vartheta}{\partial C}\f$.

\level 2

   \maintainer Amadeus Gebauer
*/

/*----------------------------------------------------------------------*/


#include "growth.H"
#include "growth_law.H"
#include "matpar_material.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/voigt_notation.H"

/*----------------------------------------------------------------------------*/
MAT::GrowthLaw::GrowthLaw() : params_(NULL) {}

/*----------------------------------------------------------------------------*/
MAT::GrowthLaw::GrowthLaw(MAT::PAR::Parameter* params) : params_(params) {}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawDyn::GrowthLawDyn(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), abstol_(matdata->GetDouble("TOL"))
{
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawDyn::CreateMaterial() { return Teuchos::null; }

/*----------------------------------------------------------------------------*/
MAT::GrowthLawDyn::GrowthLawDyn() : GrowthLaw() {}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawDyn::GrowthLawDyn(MAT::PAR::Parameter* params) : GrowthLaw(params) {}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawDyn::Evaluate(double* thetainit, const double& thetaold,
    LINALG::Matrix<6, 1>* dthetadC, MAT::Growth& matgrowth, const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  double dt = params.get<double>("delta time", -1.0);

  double theta = *thetainit;

  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatdach(true);
  LINALG::Matrix<NUM_STRESS_3D, 1> Sdachvec(true);

  // identity tensor I
  LINALG::Matrix<NUM_STRESS_3D, 1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D, 1> Cvec(*glstrain);
  Cvec.Scale(2.0);
  Cvec += Id;

  // calculate growth part F_g of the deformation gradient F
  LINALG::Matrix<3, 3> F_g(true);
  CalcFg(theta, thetaold, gp, defgrd, refdir, curdir, histdefgrd, F_g);
  // calculate F_g^(-1)
  LINALG::Matrix<3, 3> F_ginv(true);
  F_ginv.Invert(F_g);
  // elastic deformation gradient F_e = F * F_g^(-1)
  LINALG::Matrix<3, 3> defgrddach(true);   //*defgrd);
  defgrddach.MultiplyNN(*defgrd, F_ginv);  // Scale(1.0 / theta);
  // elastic right Cauchy-Green Tensor Cdach = F_e^T * F_e (= F_g^-T C F_g^-1)
  LINALG::Matrix<3, 3> Cdach(true);
  Cdach.MultiplyTN(defgrddach, defgrddach);

  // transform Cdach into a vector
  LINALG::Matrix<6, 1> Cdachvec(true);
  VStrainUtils::MatrixToVector(Cdach, Cdachvec);

  // elastic Green Lagrange strain
  LINALG::Matrix<NUM_STRESS_3D, 1> glstraindachvec(Cdachvec);
  glstraindachvec -= Id;
  glstraindachvec.Scale(0.5);
  // elastic 2 PK stress and constitutive matrix
  matgrowth.EvaluateElastic(&defgrddach, &glstraindachvec, &Sdachvec, &cmatdach, params, eleGID);

  // the growth trigger (i.e. stress, strain, ...)
  double growthtrigger = 0.0;
  EvaluateGrowthTrigger(growthtrigger, Sdachvec, Cdachvec, refdir, theta, consttrig);

  // evaluate growth function
  double growthfunc = 0.0;
  double dgrowthfunctheta = 0.0;
  EvaluateGrowthFunction(growthfunc, growthtrigger, theta);
  EvaluateGrowthFunctionDerivTheta(
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
      F_g.Scale(0.0);
      CalcFg(thetatemp, thetaold, gp, defgrd, refdir, curdir, histdefgrd, F_g);
      // calculate F_g^(-1)
      F_ginv.Scale(0.0);
      F_ginv.Invert(F_g);
      // elastic deformation gradient F_e = F * F_g^(-1)
      defgrddach.Scale(0.0);
      defgrddach.MultiplyNN(*defgrd, F_ginv);
      // elastic right Cauchy-Green tensor Cdach = F_e^T * F_e (= F_g^-T C F_g^-1)
      Cdach.Scale(0.0);
      Cdach.MultiplyTN(defgrddach, defgrddach);

      // transform Cdach into a vector
      Cdachvec.Scale(0.0);
      VStrainUtils::MatrixToVector(Cdach, Cdachvec);

      glstraindachvec = Cdachvec;
      glstraindachvec -= Id;
      glstraindachvec.Scale(0.5);
      cmatdach.Scale(0.0);
      Sdachvec.Scale(0.0);
      matgrowth.EvaluateElastic(
          &defgrddach, &glstraindachvec, &Sdachvec, &cmatdach, params, eleGID);

      // the growth trigger (i.e. stress, strain, ...)
      growthtrigger = 0.0;
      EvaluateGrowthTrigger(growthtrigger, Sdachvec, Cdachvec, refdir, thetatemp, consttrig);

      growthfunc = 0.0;
      EvaluateGrowthFunction(growthfunc, growthtrigger, thetatemp);

      residualtemp = thetaold - thetatemp + growthfunc * dt;
    }  // end of damping loop
    residual = residualtemp;
    theta = thetatemp;

    // evaluate derivative of growth function w.r.t. growth factor
    EvaluateGrowthFunctionDerivTheta(
        dgrowthfunctheta, growthtrigger, theta, Cdachvec, Sdachvec, cmatdach, F_g, refdir);

    if (omega <= omegamin and abs(residualtemp) > (1.0 - 0.5 * omega) * abs(residual))
    {
      std::cout << gp << ": Theta " << thetatemp << " residual " << residualtemp << " stress "
                << growthtrigger << std::endl;
      // dserror("no damping coefficient found");
    }

  }  // end of local Newton iteration


  if (localistep == maxstep && abs(residual) > Parameter()->abstol_)
    dserror(
        "local Newton iteration did not converge after %i steps: residual: %e, thetaold: %f,"
        " theta:  %f, growthtrigger: %e",
        maxstep, residual, thetaold, theta, growthtrigger);

  thetaquer = 1.0 - dgrowthfunctheta * dt;

  // constitutive matrix including growth cmat = F_g^-1 F_g^-1 cmatdach F_g^-T F_g^-T
  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelastic(true);
  cmatelastic = GrowthVolumetric::PullBack4Tensor(F_ginv, cmatdach);
  // calculate stress
  // 2PK stress S = F_g^-1 Sdach F_g^-T
  LINALG::Matrix<3, 3> Sdach(true);
  VStressUtils::ToMatrix(Sdachvec, Sdach);

  LINALG::Matrix<3, 3> tmp(true);
  tmp.MultiplyNT(Sdach, F_ginv);
  LINALG::Matrix<3, 3> S(true);
  S.MultiplyNN(F_ginv, tmp);

  LINALG::Matrix<6, 1> Svec(true);
  VStressUtils::MatrixToVector(S, Svec);

  LINALG::Matrix<NUM_STRESS_3D, 1> dgrowthfuncdCvec(true);
  EvaluateGrowthFunctionDerivC(
      dgrowthfuncdCvec, growthtrigger, theta, Cvec, Svec, cmatelastic, refdir);

  // d theta/ d C
  dthetadC->Update(dt / thetaquer, dgrowthfuncdCvec);

  *thetainit = theta;

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawDyn::EvaluatePDeriv(double* thetainit, const double& thetaold,
    Teuchos::RCP<MAT::So3Material> matelastic, const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params, const int eleGID)
{
  dserror("This must still be implemented for dynamic type growth laws.");
}

/*----------------------------------------------------------------------------*/
double MAT::GrowthLawDyn::DensityScale(const double theta)
{
  // nothing to do
  return 0;
}
double MAT::GrowthLawDyn::DensityDerivScale(const double theta)
{
  // nothing to do
  return 0;
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawStatic::GrowthLawStatic() : GrowthLaw() {}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawStatic::GrowthLawStatic(MAT::PAR::Parameter* params) : GrowthLaw(params) {}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawStatic::EvaluatePDeriv(double* thetainit, const double& thetaold,
    Teuchos::RCP<MAT::So3Material> matelastic, const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params, const int eleGID)
{
  dserror("Use implementation of a specific growth law overriding this function.");
}

/*----------------------------------------------------------------------------*/
double MAT::GrowthLawStatic::DensityScale(const double theta)
{
  // nothing to do
  return 0;
}
/*----------------------------------------------------------------------------*/
double MAT::GrowthLawStatic::DensityDerivScale(const double theta)
{
  // nothing to do
  return 0;
}



/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawAnisoStrain::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawAnisoStrain::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawAnisoStrain(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawAnisoStrain::GrowthLawAnisoStrain(Teuchos::RCP<MAT::PAR::Material> matdata)
    : GrowthLawDyn(matdata),
      tau_(matdata->GetDouble("TAU")),
      taurev_(matdata->GetDouble("TAU_REV")),
      thetamin_(matdata->GetDouble("THETA_MIN")),
      thetamax_(matdata->GetDouble("THETA_MAX")),
      gamma_(matdata->GetDouble("GAMMA")),
      gammarev_(matdata->GetDouble("GAMMA_REV")),
      lambdacrit_(matdata->GetDouble("LAMBDA_CRIT"))
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawAnisoStrain::GrowthLawAnisoStrain() : GrowthLawDyn(NULL) {}


/*----------------------------------------------------------------------------*/
MAT::GrowthLawAnisoStrain::GrowthLawAnisoStrain(MAT::PAR::GrowthLawAnisoStrain* params)
    : GrowthLawDyn(params)
{
}



/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStrain::EvaluateGrowthFunction(
    double& growthfunc, const double growthtrig, const double theta)
{
  MAT::PAR::GrowthLawAnisoStrain* params = Parameter();
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

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStrain::EvaluateGrowthTrigger(double& growthtrig,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec, const LINALG::Matrix<6, 1>& Cdachvec,
    const LINALG::Matrix<3, 1>& direction, const double& theta, const double& consttrig)
{
  // transform Cdachvec into a matrix
  LINALG::Matrix<3, 3> Cdach(true);
  VStrainUtils::ToMatrix(Cdachvec, Cdach);

  LINALG::Matrix<3, 1> CdachDir(true);
  CdachDir.MultiplyNN(1.0, Cdach, direction);

  // elastic fiber stretch
  growthtrig =
      sqrt(CdachDir(0) * direction(0) + CdachDir(1) * direction(1) + CdachDir(2) * direction(2));

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStrain::EvaluateGrowthFunctionDerivTheta(double& dgrowthfunctheta,
    double growthtrig, double theta, const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const LINALG::Matrix<3, 3>& F_g, const LINALG::Matrix<3, 1>& direction)
{
  MAT::PAR::GrowthLawAnisoStrain* params = Parameter();
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

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStrain::EvaluateGrowthFunctionDerivC(
    LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double growthtrig, double theta,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec, const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat, const LINALG::Matrix<3, 1>& direction)
{
  MAT::PAR::GrowthLawAnisoStrain* params = Parameter();
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


  LINALG::Matrix<3, 3> dgrowthfuncdC(true);
  dgrowthfuncdC.MultiplyNT(direction, direction);

  dgrowthfuncdC.Scale(ktheta / (theta * theta * 2. * growthtrig));

  // transform dgrowthfuncdC into a vector in Voigt notation
  dgrowthfuncdCvec(0) = dgrowthfuncdC(0, 0);
  dgrowthfuncdCvec(1) = dgrowthfuncdC(1, 1);
  dgrowthfuncdCvec(2) = dgrowthfuncdC(2, 2);
  dgrowthfuncdCvec(3) = dgrowthfuncdC(0, 1);
  dgrowthfuncdCvec(4) = dgrowthfuncdC(1, 2);
  dgrowthfuncdCvec(5) = dgrowthfuncdC(0, 2);

  return;
}

///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStrain::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, LINALG::Matrix<3, 3>& F_g)
{
  LINALG::Matrix<3, 3> fdf(true);
  fdf.MultiplyNT(refdir, refdir);

  for (int i = 0; i < 3; i++) F_g(i, i) = 1.0;
  F_g.Update(theta - 1.0, fdf, 1.0);

  return;
}


/*----------------------------------------------------------------------------*/
double MAT::GrowthLawAnisoStrain::DensityScale(const double theta)
{
  // in anisotropic uni-directional growth case: rho_0 = theta * rho_dach,
  // thus d(rho_0)/d(theta) = rho_dach
  return theta;
}
/*----------------------------------------------------------------------------*/
double MAT::GrowthLawAnisoStrain::DensityDerivScale(const double theta) { return 1.0; }


/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawAnisoStress::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawAnisoStress::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawAnisoStress(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawAnisoStress::GrowthLawAnisoStress(Teuchos::RCP<MAT::PAR::Material> matdata)
    : GrowthLawDyn(matdata),
      tau_(matdata->GetDouble("TAU")),
      taurev_(matdata->GetDouble("TAU_REV")),
      thetamin_(matdata->GetDouble("THETA_MIN")),
      thetamax_(matdata->GetDouble("THETA_MAX")),
      gamma_(matdata->GetDouble("GAMMA")),
      gammarev_(matdata->GetDouble("GAMMA_REV")),
      pcrit_(matdata->GetDouble("P_CRIT"))
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawAnisoStress::GrowthLawAnisoStress() : GrowthLawDyn(NULL) {}


/*----------------------------------------------------------------------------*/
MAT::GrowthLawAnisoStress::GrowthLawAnisoStress(MAT::PAR::GrowthLawAnisoStress* params)
    : GrowthLawDyn(params)
{
}



/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStress::EvaluateGrowthFunction(
    double& growthfunc, const double growthtrig, const double theta)
{
  MAT::PAR::GrowthLawAnisoStress* params = Parameter();
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

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStress::EvaluateGrowthTrigger(double& growthtrig,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec, const LINALG::Matrix<6, 1>& Cdachvec,
    const LINALG::Matrix<3, 1>& direction, const double& theta, const double& consttrig)
{
  // trace of elastic Mandel stress tensor
  growthtrig = Cdachvec(0) * Sdachvec(0) + Cdachvec(1) * Sdachvec(1) + Cdachvec(2) * Sdachvec(2) +
               Cdachvec(3) * Sdachvec(3) + Cdachvec(4) * Sdachvec(4) + Cdachvec(5) * Sdachvec(5);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStress::EvaluateGrowthFunctionDerivTheta(double& dgrowthfunctheta,
    double growthtrig, double theta, const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const LINALG::Matrix<3, 3>& F_g, const LINALG::Matrix<3, 1>& direction)
{
  MAT::PAR::GrowthLawAnisoStress* params = Parameter();
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
    ;
  }

  // calculate F_g^(-1)
  LINALG::Matrix<3, 3> F_ginv(true);
  F_ginv.Invert(F_g);

  // calculate d(F_g)/d(theta) = 1 - f_0 \otimes f_0
  LINALG::Matrix<3, 3> fdf(true);
  fdf.MultiplyNT(direction, direction);

  LINALG::Matrix<3, 3> dFgdtheta(true);
  for (int i = 0; i < 3; i++) dFgdtheta(i, i) = 1.0;

  dFgdtheta.Update(-1.0, fdf, 1.0);


  // transform Cdachvec into a matrix
  LINALG::Matrix<3, 3> Cdach(true);
  VStrainUtils::ToMatrix(Cdachvec, Cdach);

  LINALG::Matrix<3, 3> dFgT_Cdach(true);
  dFgT_Cdach.MultiplyTN(dFgdtheta, Cdach);
  LINALG::Matrix<3, 3> dCdachdtheta(true);
  dCdachdtheta.MultiplyTN(F_ginv, dFgT_Cdach);

  LINALG::Matrix<3, 3> dFg_Fginv(true);
  dFg_Fginv.MultiplyNN(dFgdtheta, F_ginv);
  LINALG::Matrix<3, 3> Cdach_dFg_Fginv(true);
  Cdach_dFg_Fginv.MultiplyNN(Cdach, dFg_Fginv);

  dCdachdtheta.Update(-1.0, Cdach_dFg_Fginv, -1.0);

  // transform dCdachdtheta into a vector
  LINALG::Matrix<6, 1> dCdachdthetavec(true);
  VStrainUtils::MatrixToVector(dCdachdtheta, dCdachdthetavec);

  LINALG::Matrix<6, 1> dSdachdthetavec(true);

  for (int i = 0; i < 6; i++)
  {
    dSdachdthetavec(i) =
        0.5 * (cmatelastic(i, 0) * dCdachdthetavec(0) + cmatelastic(i, 1) * dCdachdthetavec(1) +
                  cmatelastic(i, 2) * dCdachdthetavec(2) + cmatelastic(i, 3) * dCdachdthetavec(3) +
                  cmatelastic(i, 4) * dCdachdthetavec(4) + cmatelastic(i, 5) * dCdachdthetavec(5));
  }

  // derivative of trace of Mandel stress tensor w.r.t. theta
  // cf. Göktepe et al. J Theor Biol (2010), eq. 26 ff.
  const double dgrowthtrig = dCdachdthetavec(0) * Sdachvec(0) + dCdachdthetavec(1) * Sdachvec(1) +
                             dCdachdthetavec(2) * Sdachvec(2) + dCdachdthetavec(3) * Sdachvec(3) +
                             dCdachdthetavec(4) * Sdachvec(4) + dCdachdthetavec(5) * Sdachvec(5)

                             + dSdachdthetavec(0) * Cdachvec(0) + dSdachdthetavec(1) * Cdachvec(1) +
                             dSdachdthetavec(2) * Cdachvec(2) + dSdachdthetavec(3) * Cdachvec(3) +
                             dSdachdthetavec(4) * Cdachvec(4) + dSdachdthetavec(5) * Cdachvec(5);


  dgrowthfunctheta = dktheta * (growthtrig - pcrit) + ktheta * dgrowthtrig;

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStress::EvaluateGrowthFunctionDerivC(
    LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double traceM, double theta,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec, const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat, const LINALG::Matrix<3, 1>& direction)
{
  MAT::PAR::GrowthLawAnisoStress* params = Parameter();
  // parameters
  const double tau = params->tau_;
  const double taurev = params->taurev_;
  const double thetamin = params->thetamin_;
  const double thetamax = params->thetamax_;
  const double gamma = params->gamma_;
  const double gammarev = params->gammarev_;
  const double pcrit = params->pcrit_;

  double ktheta = 0.;

  if (traceM >= pcrit)
    ktheta = (1. / tau) * pow(((thetamax - theta) / (thetamax - thetamin)), gamma);
  else
    ktheta = (1. / taurev) * pow(((theta - thetamin) / (thetamax - thetamin)), gammarev);


  for (int j = 0; j < NUM_STRESS_3D; j++)
  {
    double Ccmatelasj = cmat(0, j) * Cvec(0) + cmat(1, j) * Cvec(1) + cmat(2, j) * Cvec(2) +
                        cmat(3, j) * Cvec(3) + cmat(4, j) * Cvec(4) + cmat(5, j) * Cvec(5);

    dgrowthfuncdCvec(j) = ktheta * (Svec(j) + 0.5 * Ccmatelasj);
  }

  return;
}


///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStress::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, LINALG::Matrix<3, 3>& F_g)
{
  LINALG::Matrix<3, 3> fdf(true);
  fdf.MultiplyNT(refdir, refdir);

  for (int i = 0; i < 3; i++) F_g(i, i) = theta;
  F_g.Update(1.0 - theta, fdf, 1.0);

  return;
}


/*----------------------------------------------------------------------------*/
double MAT::GrowthLawAnisoStress::DensityScale(const double theta)
{
  // in anisotropic bi-directional growth case: rho_0 = theta^2 * rho_dach,
  // thus d(rho_0)/d(theta) = 2 * theta * rho_dach
  return theta * theta;
}
/*----------------------------------------------------------------------------*/
double MAT::GrowthLawAnisoStress::DensityDerivScale(const double theta) { return 2.0 * theta; }



/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawAnisoStrainConstTrig::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawAnisoStrainConstTrig::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawAnisoStrainConstTrig(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawAnisoStrainConstTrig::GrowthLawAnisoStrainConstTrig(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : GrowthLawAnisoStrain(matdata)
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawAnisoStrainConstTrig::GrowthLawAnisoStrainConstTrig() : GrowthLawAnisoStrain(NULL) {}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawAnisoStrainConstTrig::GrowthLawAnisoStrainConstTrig(
    MAT::PAR::GrowthLawAnisoStrainConstTrig* params)
    : GrowthLawAnisoStrain(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStrainConstTrig::EvaluateGrowthTrigger(double& growthtrig,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec, const LINALG::Matrix<6, 1>& Cdachvec,
    const LINALG::Matrix<3, 1>& direction, const double& theta, const double& consttrig)
{
  growthtrig = consttrig;

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStrainConstTrig::EvaluateGrowthFunctionDerivTheta(double& dgrowthfunctheta,
    double growthtrig, double theta, const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const LINALG::Matrix<3, 3>& F_g, const LINALG::Matrix<3, 1>& direction)
{
  MAT::PAR::GrowthLawAnisoStrainConstTrig* params = Parameter();
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

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStrainConstTrig::EvaluateGrowthFunctionDerivC(
    LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double traceM, double theta,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec, const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat, const LINALG::Matrix<3, 1>& direction)
{
  for (int j = 0; j < NUM_STRESS_3D; j++)
  {
    dgrowthfuncdCvec(j) = 0.0;
  }

  return;
}



/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawAnisoStressConstTrig::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawAnisoStressConstTrig::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawAnisoStressConstTrig(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawAnisoStressConstTrig::GrowthLawAnisoStressConstTrig(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : GrowthLawAnisoStress(matdata)
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawAnisoStressConstTrig::GrowthLawAnisoStressConstTrig() : GrowthLawAnisoStress(NULL) {}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawAnisoStressConstTrig::GrowthLawAnisoStressConstTrig(
    MAT::PAR::GrowthLawAnisoStressConstTrig* params)
    : GrowthLawAnisoStress(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStressConstTrig::EvaluateGrowthTrigger(double& growthtrig,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec, const LINALG::Matrix<6, 1>& Cdachvec,
    const LINALG::Matrix<3, 1>& direction, const double& theta, const double& consttrig)
{
  growthtrig = consttrig;

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStressConstTrig::EvaluateGrowthFunctionDerivTheta(double& dgrowthfunctheta,
    double growthtrig, double theta, const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const LINALG::Matrix<3, 3>& F_g, const LINALG::Matrix<3, 1>& direction)
{
  MAT::PAR::GrowthLawAnisoStressConstTrig* params = Parameter();
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
    ;
  }

  dgrowthfunctheta = dktheta * (growthtrig - pcrit);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAnisoStressConstTrig::EvaluateGrowthFunctionDerivC(
    LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double traceM, double theta,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec, const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat, const LINALG::Matrix<3, 1>& direction)
{
  for (int j = 0; j < NUM_STRESS_3D; j++)
  {
    dgrowthfuncdCvec(j) = 0.0;
  }

  return;
}



/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawIsoStress::CreateMaterial() { return Teuchos::null; }

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawIsoStress::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawIsoStress(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawIsoStress::GrowthLawIsoStress(Teuchos::RCP<MAT::PAR::Material> matdata)
    : GrowthLawDyn(matdata),
      thetaplus_(matdata->GetDouble("THETAPLUS")),
      kthetaplus_(matdata->GetDouble("KPLUS")),
      mthetaplus_(matdata->GetDouble("MPLUS")),
      thetaminus_(matdata->GetDouble("THETAMINUS")),
      kthetaminus_(matdata->GetDouble("KMINUS")),
      mthetaminus_(matdata->GetDouble("MMINUS")),
      hommandel_(matdata->GetDouble("HOMMANDEL"))
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawIsoStress::GrowthLawIsoStress() : GrowthLawDyn(NULL) {}


/*----------------------------------------------------------------------------*/
MAT::GrowthLawIsoStress::GrowthLawIsoStress(MAT::PAR::GrowthLawIsoStress* params)
    : GrowthLawDyn(params)
{
}



/*----------------------------------------------------------------------------*/
void MAT::GrowthLawIsoStress::EvaluateGrowthFunction(
    double& growthfunc, const double traceM, const double theta)
{
  MAT::PAR::GrowthLawIsoStress* params = Parameter();
  // parameters
  const double hommandel = params->hommandel_;
  const double thetaplus = params->thetaplus_;      // 1.5;
  const double kthetaplus = params->kthetaplus_;    // 1.0; 0.5;
  const double mthetaplus = params->mthetaplus_;    // 2.0; 4.0;
  const double thetaminus = params->thetaminus_;    // 0.5;
  const double kthetaminus = params->kthetaminus_;  // 2.0; 0.25;
  const double mthetaminus = params->mthetaminus_;  // 3.0; 5.0;

  double ktheta = 0.0;

  if (traceM >= hommandel)
    ktheta = kthetaplus * pow((thetaplus - theta) / (thetaplus - 1.0), mthetaplus);
  else
    ktheta = kthetaminus * pow((theta - thetaminus) / (1.0 - thetaminus), mthetaminus);


  growthfunc = ktheta * (traceM - hommandel);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawIsoStress::EvaluateGrowthTrigger(double& growthtrig,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec, const LINALG::Matrix<6, 1>& Cdachvec,
    const LINALG::Matrix<3, 1>& direction, const double& theta, const double& consttrig)
{
  // trace of elastic Mandel stress tensor
  growthtrig = Cdachvec(0) * Sdachvec(0) + Cdachvec(1) * Sdachvec(1) + Cdachvec(2) * Sdachvec(2) +
               Cdachvec(3) * Sdachvec(3) + Cdachvec(4) * Sdachvec(4) + Cdachvec(5) * Sdachvec(5);

  return;
}


/*----------------------------------------------------------------------------*/
void MAT::GrowthLawIsoStress::EvaluateGrowthFunctionDerivTheta(double& dgrowthfunctheta,
    double traceM, double theta, const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Sdachvec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic,
    const LINALG::Matrix<3, 3>& F_g, const LINALG::Matrix<3, 1>& direction)
{
  MAT::PAR::GrowthLawIsoStress* params = Parameter();
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

  if (traceM >= hommandel)
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

  const double dtraceM = -(2.0 * traceM + temp) / theta;

  dgrowthfunctheta = (dktheta * (traceM - hommandel) + ktheta * dtraceM);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawIsoStress::EvaluateGrowthFunctionDerivC(
    LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdCvec, double growthtrig, double theta,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Cvec, const LINALG::Matrix<NUM_STRESS_3D, 1>& Svec,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat, const LINALG::Matrix<3, 1>& direction)
{
  MAT::PAR::GrowthLawIsoStress* params = Parameter();
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

  return;
}


///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void MAT::GrowthLawIsoStress::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, LINALG::Matrix<3, 3>& F_g)
{
  for (int i = 0; i < 3; i++) F_g(i, i) = theta;

  return;
}


/*----------------------------------------------------------------------------*/
double MAT::GrowthLawIsoStress::DensityScale(const double theta)
{
  // in isotropic growth case: rho_0 = theta^3 * rho_dach,
  // thus d(rho_0)/d(theta) = 3 * theta^2 * rho_dach
  return theta * theta * theta;
}
/*----------------------------------------------------------------------------*/
double MAT::GrowthLawIsoStress::DensityDerivScale(const double theta)
{
  return 3.0 * theta * theta;
}



/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawAC::GrowthLawAC(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      Sc1_(matdata->GetInt("SCALAR1")),
      alpha_(matdata->GetDouble("ALPHA")),
      Sc2_(matdata->GetInt("SCALAR2")),
      beta_(matdata->GetDouble("BETA"))
{
  if (Sc1_ < 1) dserror("At least on scalar field must induce growth");
  if (alpha_ < 0) dserror("The influence of scalar field SCALAR1 to growth can't be negativ");
  if (Sc2_ < 1)
    dserror("If you choose a second scalar field to induce growth, choose a existing one!!");
  if (beta_ < 0) dserror("The influence of scalar field SCALAR2 to growth can't be negativ");
}

Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawAC::CreateMaterial() { return Teuchos::null; }

Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawAC::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawAC(this));
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawAC::GrowthLawAC() : GrowthLawStatic(NULL) {}


/*----------------------------------------------------------------------------*/
MAT::GrowthLawAC::GrowthLawAC(MAT::PAR::GrowthLawAC* params) : GrowthLawStatic(params) {}



/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAC::Evaluate(double* theta, const double& thetaold,
    LINALG::Matrix<6, 1>* dthetadC, MAT::Growth& matgrowth, const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int eleGID)
{
  // get Gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("No Gauss point number provided in material.");
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

  LINALG::Matrix<3, 3> C(true);
  C.MultiplyTN(*defgrd, *defgrd);
  LINALG::Matrix<3, 3> Cinv(true);
  Cinv.Invert(C);

  // linearization of growth law
  LINALG::Matrix<3, 3> dThetadC(Cinv);
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
  Teuchos::RCP<LINALG::Matrix<6, 1>> dlambda_dC = Teuchos::rcp(new LINALG::Matrix<6, 1>(true));
  (*dlambda_dC)(0) = Cinv(0, 0);
  (*dlambda_dC)(1) = Cinv(1, 1);
  (*dlambda_dC)(2) = Cinv(2, 2);
  (*dlambda_dC)(3) = Cinv(0, 1);
  (*dlambda_dC)(4) = Cinv(1, 2);
  (*dlambda_dC)(5) = Cinv(0, 2);

  // save values in parameter list
  params.set<double>("lambda", lambda);
  params.set<Teuchos::RCP<LINALG::Matrix<6, 1>>>("dlambda_dC", dlambda_dC);
}


///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void MAT::GrowthLawAC::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, LINALG::Matrix<3, 3>& F_g)
{
  for (int i = 0; i < 3; i++) F_g(i, i) = theta;

  return;
}


/*----------------------------------------------------------------------------*/
double MAT::GrowthLawAC::DensityScale(const double theta)
{
  // in isotropic growth case: rho_0 = theta^3 * rho_dach,
  // thus d(rho_0)/d(theta) = 3 * theta^2 * rho_dach
  return theta * theta * theta;
}
/*----------------------------------------------------------------------------*/
double MAT::GrowthLawAC::DensityDerivScale(const double theta) { return 3.0 * theta * theta; }



/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawACRadial::GrowthLawACRadial(Teuchos::RCP<MAT::PAR::Material> matdata)
    : GrowthLawAC(matdata)
{
}

Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawACRadial::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawACRadial(this));
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawACRadial::GrowthLawACRadial() : GrowthLawStatic(NULL) {}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawACRadial::GrowthLawACRadial(MAT::PAR::GrowthLawAC* params) : GrowthLawStatic(params)
{
}



/*----------------------------------------------------------------------------*/
void MAT::GrowthLawACRadial::Evaluate(double* theta, const double& thetaold,
    LINALG::Matrix<6, 1>* dthetadC, MAT::Growth& matgrowth, const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int eleGID)
{
  // get Gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("No Gauss point number provided in material.");
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

  LINALG::Matrix<3, 3> C(true);
  C.MultiplyTN(*defgrd, *defgrd);
  LINALG::Matrix<3, 3> Cinv(true);
  Cinv.Invert(C);

  // linearization of growth law
  LINALG::Matrix<3, 3> dThetadC(Cinv);
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
  Teuchos::RCP<LINALG::Matrix<6, 1>> dlambda_dC = Teuchos::rcp(new LINALG::Matrix<6, 1>(true));
  (*dlambda_dC)(0) = Cinv(0, 0);
  (*dlambda_dC)(1) = Cinv(1, 1);
  (*dlambda_dC)(2) = Cinv(2, 2);
  (*dlambda_dC)(3) = Cinv(0, 1);
  (*dlambda_dC)(4) = Cinv(1, 2);
  (*dlambda_dC)(5) = Cinv(0, 2);

  // save values in parameter list
  params.set<double>("lambda", lambda);
  params.set<Teuchos::RCP<LINALG::Matrix<6, 1>>>("dlambda_dC", dlambda_dC);
}



///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void MAT::GrowthLawACRadial::CalcFg(
    const double& theta, const double& thetaold, const int& gp, const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<3, 1>& refdir, const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, LINALG::Matrix<3, 3>& F_g

)
{
  LINALG::Matrix<3, 3> ndn(true);
  ndn.MultiplyNT(curdir.at(gp), curdir.at(gp));

  LINALG::Matrix<3, 3> F_g_incr(true);
  for (int i = 0; i < 3; i++) F_g_incr(i, i) = 1.0;
  F_g_incr.Update((theta - thetaold) / thetaold, ndn, 1.0);

  F_g.MultiplyNN(F_g_incr, histdefgrd.at(gp));

  return;
}



/*----------------------------------------------------------------------------*/
double MAT::GrowthLawACRadial::DensityScale(const double theta)
{
  // in anisotropic uni-directional growth case: rho_0 = theta * rho_dach,
  // thus d(rho_0)/d(theta) = rho_dach
  return theta;
}
/*----------------------------------------------------------------------------*/
double MAT::GrowthLawACRadial::DensityDerivScale(const double theta) { return 1.0; }



/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawACRadialRefConc::GrowthLawACRadialRefConc(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : GrowthLawAC(matdata)
{
}

Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawACRadialRefConc::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawACRadialRefConc(this));
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawACRadialRefConc::GrowthLawACRadialRefConc() : GrowthLawStatic(NULL) {}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawACRadialRefConc::GrowthLawACRadialRefConc(MAT::PAR::GrowthLawAC* params)
    : GrowthLawStatic(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawACRadialRefConc::Evaluate(double* theta, const double& thetaold,
    LINALG::Matrix<6, 1>* dthetadC, MAT::Growth& matgrowth, const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  // get pointer to vector containing the scalar values at the Gauß points
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
  //  params.set< Teuchos::RCP<LINALG::Matrix<6,1> > >("dlambda_dC",dlambda_dC);
}


///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void MAT::GrowthLawACRadialRefConc::CalcFg(const double& theta, const double& thetaold,
    const int& gp, const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, LINALG::Matrix<3, 3>& F_g)
{
  LINALG::Matrix<3, 3> ndn(true);
  ndn.MultiplyNT(curdir.at(gp), curdir.at(gp));

  LINALG::Matrix<3, 3> F_g_incr(true);
  for (int i = 0; i < 3; i++) F_g_incr(i, i) = 1.0;
  F_g_incr.Update((theta - thetaold) / thetaold, ndn, 1.0);

  F_g.MultiplyNN(F_g_incr, histdefgrd.at(gp));

  return;
}


/*----------------------------------------------------------------------------*/
double MAT::GrowthLawACRadialRefConc::DensityScale(const double theta)
{
  // in anisotropic uni-directional growth case: rho_0 = theta * rho_dach,
  // thus d(rho_0)/d(theta) = rho_dach
  return theta;
}
/*----------------------------------------------------------------------------*/
double MAT::GrowthLawACRadialRefConc::DensityDerivScale(const double theta) { return 1.0; }



/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawConst::CreateMaterial() { return Teuchos::null; }

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawConst::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawConst(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawConst::GrowthLawConst(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  Epetra_Map dummy_map(1, 1, 0, *(DRT::Problem::Instance()->GetNPGroup()->LocalComm()));
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(thetarate)->PutScalar(matdata->GetDouble("THETARATE"));
}

void MAT::PAR::GrowthLawConst::OptParams(std::map<std::string, int>* pnames)
{
  pnames->insert(std::pair<std::string, int>("THETARATE", thetarate));
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawConst::GrowthLawConst() : GrowthLawStatic(NULL) {}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawConst::GrowthLawConst(MAT::PAR::GrowthLawConst* params) : GrowthLawStatic(params) {}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawConst::Evaluate(double* theta, const double& thetaold,
    LINALG::Matrix<6, 1>* dthetadC, MAT::Growth& matgrowth, const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, const double& consttrig,
    Teuchos::ParameterList& params, const int eleGID)
{
  double dt = params.get<double>("delta time", -1.0);
  // map in GetParameter can now calculate LID, so we do not need it here       05/2017 birzle
  // int eleID = DRT::Problem::Instance()->GetDis("structure")->ElementColMap()->LID(eleGID);
  *theta = thetaold + Parameter()->GetParameter(Parameter()->thetarate, eleGID) * dt;
  //*theta = Parameter()->GetParameter(Parameter()->thetarate,eleID);

  dthetadC->Scale(0.0);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawConst::EvaluatePDeriv(double* theta, const double& thetaold,
    Teuchos::RCP<MAT::So3Material> matelastic, const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params, const int eleGID)
{
  double dt = params.get<double>("delta time", -1.0);

  int deriv = params.get<int>("matparderiv", -1);
  if (deriv == Parameter()->thetarate)
  {
    *theta = dt;
  }
  else if (deriv == -1)
    dserror("you should only end up here with a valid matparderiv flag!");
}

///*----------------------------------------------------------------------*
// | calculate growth part of deformation gradient                mh 04/17|
// *----------------------------------------------------------------------*/
void MAT::GrowthLawConst::CalcFg(const double& theta, const double& thetaold, const int& gp,
    const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<3, 1>& refdir,
    const std::vector<LINALG::Matrix<3, 1>>& curdir,
    const std::vector<LINALG::Matrix<3, 3>>& histdefgrd, LINALG::Matrix<3, 3>& F_g)
{
  for (int i = 0; i < 3; i++) F_g(i, i) = theta;

  return;
}


/*----------------------------------------------------------------------------*/
double MAT::GrowthLawConst::DensityScale(const double theta)
{
  // in isotropic growth case: rho_0 = theta^3 * rho_dach,
  // thus d(rho_0)/d(theta) = 3 * theta^2 * rho_dach
  return theta * theta * theta;
}
/*----------------------------------------------------------------------------*/
double MAT::GrowthLawConst::DensityDerivScale(const double theta) { return 3.0 * theta * theta; }
