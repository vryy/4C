/*----------------------------------------------------------------------*/
/*!
 * \file growth_law.cpp
 *
 * <pre>
 * Maintainer: Moritz Thon
 *             thon@lnm.mw.tum.de
 *             http://www.mhpc.mw.tum.de
 *             089 - 289-10364
 * </pre>
 *----------------------------------------------------------------------*/


#include "growth_law.H"
#include "matpar_material.H"


/*----------------------------------------------------------------------------*/
MAT::GrowthLaw::GrowthLaw()
  : params_(NULL),
    cfac_(1,1.0)
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLaw::GrowthLaw(MAT::PAR::Parameter* params)
  : params_(params),
    cfac_(1,1.0)
{
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawDyn::GrowthLawDyn(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  abstol_(matdata->GetDouble("TOL"))
{
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawDyn::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawDyn::GrowthLawDyn()
  : GrowthLaw()
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawDyn::GrowthLawDyn(MAT::PAR::Parameter* params)
  : GrowthLaw(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawDyn::Evaluate(double* thetainit,
                                 const double& thetaold,
                                 LINALG::Matrix<6,1>* dthetadC,
                                 Teuchos::RCP<MAT::So3Material> matelastic,
                                 const LINALG::Matrix<3,3>* defgrd,
                                 const LINALG::Matrix<6,1>* glstrain,
                                 Teuchos::ParameterList& params,
                                 const int eleGID)
{

  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1)
    dserror("no Gauss point number provided in material");

  double dt = params.get<double>("delta time", -1.0);

  double theta = *thetainit;

  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelastic(true);
  LINALG::Matrix<NUM_STRESS_3D, 1> Sdach(true);

  // identity tensor I
  LINALG::Matrix<NUM_STRESS_3D, 1> Id(true);
  for (int i = 0; i < 3; i++)
    Id(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D, 1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  // elastic right Cauchy-Green Tensor Cdach = F_g^-T C F_g^-1
  LINALG::Matrix<NUM_STRESS_3D, 1> Cdach(C);
  Cdach.Scale(1.0 / theta / theta);
  LINALG::Matrix<3, 3> defgrddach(*defgrd);
  defgrddach.Scale(1.0 / theta);
  // elastic Green Lagrange strain
  LINALG::Matrix<NUM_STRESS_3D, 1> glstraindach(Cdach);
  glstraindach -= Id;
  glstraindach.Scale(0.5);
  // elastic 2 PK stress and constitutive matrix
  matelastic->Evaluate(&defgrddach,
      &glstraindach,
      params,
      &Sdach,
      &cmatelastic,
      eleGID);

  // trace of elastic Mandel stress Mdach = Cdach Sdach
  double mandel =   Cdach(0) * Sdach(0)
                          + Cdach(1) * Sdach(1)
                          + Cdach(2) * Sdach(2)
                          + Cdach(3) * Sdach(3)
                          + Cdach(4) * Sdach(4)
                          + Cdach(5) * Sdach(5);

  // evaluate growth function
  double growthfunc = 0.0;
  double dgrowthfunctheta = 0.0;
  EvaluateGrowthFunction(growthfunc, mandel, theta);
  EvaluateGrowthFunctionDerivTheta(dgrowthfunctheta, mandel, theta, Cdach, cmatelastic);

  double residual = thetaold - theta + growthfunc * dt;

  int localistep = 0;
  double thetaquer = 0.0;
  int maxstep = 30;
  // TODO: put newton tol in growthlaw material line definition

  // local Newton iteration to obtain exact theta
  while (abs(residual) > Parameter()->abstol_ && localistep < maxstep)
  {
    localistep += 1;

    //evaluate derivative of growth function w.r.t. growth factor
    thetaquer =   1.0 - dgrowthfunctheta * dt;

    // damping strategy
    double omega = 2.0;
    double thetatemp = theta;
    double residualtemp = residual;
    double omegamin = 1.0 / 64.0;
    while ( abs(residualtemp) > (1.0 - 0.5 * omega) * abs(residual) and
        omega > omegamin )
    {
      // update of theta
      omega = 0.5*omega;
      thetatemp = theta + omega * residual / thetaquer;

      // update elastic variables
      Cdach = C;
      Cdach.Scale(1.0 / thetatemp / thetatemp);
      LINALG::Matrix<3, 3> defgrddach(*defgrd);
      defgrddach.Scale(1.0 / thetatemp);
      glstraindach = Cdach;
      glstraindach -= Id;
      glstraindach.Scale(0.5);
      cmatelastic.Scale(0.0);
      Sdach.Scale(0.0);
      matelastic->Evaluate(&defgrddach,
          &glstraindach,
          params,
          &Sdach,
          &cmatelastic,
          eleGID);

      // trace of mandel stress
      mandel =   Cdach(0) * Sdach(0) + Cdach(1) * Sdach(1) + Cdach(2) * Sdach(2)
      + Cdach(3) * Sdach(3) + Cdach(4) * Sdach(4) + Cdach(5) * Sdach(5);
      //if (signmandel*mandel < 0) signmandel = -1.0*signmandel;

      growthfunc = 0.0;
      EvaluateGrowthFunction(growthfunc, mandel, thetatemp);

      residualtemp =    thetaold - thetatemp + growthfunc * dt;
    } // end of damping loop
    residual = residualtemp;
    theta = thetatemp;

    //evaluate derivative of growth function w.r.t. growth factor
    EvaluateGrowthFunctionDerivTheta(dgrowthfunctheta, mandel, theta, Cdach, cmatelastic);

    if ( omega <= omegamin and
        abs(residualtemp) > (1.0 - 0.5 * omega) * abs(residual))
    {
      std::cout << gp << ": Theta " << thetatemp << " residual "
          << residualtemp << " stress " << mandel << std::endl;
      //dserror("no damping coefficient found");
    }

  } // end of local Newton iteration

  if (localistep == maxstep && abs(residual) > Parameter()->abstol_)
    dserror("local Newton iteration did not converge after %i steps: residual: %e, thetaold: %f,"
        " theta:  %f, mandel: %e",maxstep, residual, thetaold, theta, mandel);

  thetaquer = 1.0 - dgrowthfunctheta * dt;

  // constitutive matrix including growth
  cmatelastic.Scale(1.0 / theta / theta / theta / theta);

  // 2PK stress S = F_g^-1 Sdach F_g^-T
  LINALG::Matrix<NUM_STRESS_3D, 1> S(Sdach);
  S.Scale(1.0 / theta / theta);

  LINALG::Matrix<NUM_STRESS_3D, 1> dgrowthfuncdC(true);
  EvaluateGrowthFunctionDerivC(dgrowthfuncdC,mandel,theta,C,S,cmatelastic);

  // d theta/ d C
  dthetadC->Update(dt/thetaquer,dgrowthfuncdC);

  *thetainit=theta;

}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawDyn::EvaluateNonLinMass(double* theta,
                                           const double& thetaold,
                                           LINALG::Matrix<6,1>* linmass_disp,
                                           Teuchos::RCP<MAT::So3Material> matelastic,
                                           const LINALG::Matrix<3,3>* defgrd,
                                           const LINALG::Matrix<6,1>* glstrain,
                                           Teuchos::ParameterList& params,
                                           const int eleGID )
{
  Evaluate(theta,thetaold,linmass_disp,matelastic,defgrd,glstrain,params,eleGID);
  linmass_disp->Scale(3.0*(*theta)*(*theta)*matelastic->Density());
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawStatic::GrowthLawStatic()
  : GrowthLaw()
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawStatic::GrowthLawStatic(MAT::PAR::Parameter* params)
  : GrowthLaw(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawStatic::EvaluateNonLinMass(double* theta,
                                              const double& thetaold,
                                              LINALG::Matrix<6,1>* linmass_disp,
                                              Teuchos::RCP<MAT::So3Material> matelastic,
                                              const LINALG::Matrix<3,3>* defgrd,
                                              const LINALG::Matrix<6,1>* glstrain,
                                              Teuchos::ParameterList& params,
                                              const int eleGID )
{
  // this is a static growth law
  linmass_disp->Clear();
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawLinear::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawLinear::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawLinear(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawLinear::GrowthLawLinear(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
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
MAT::GrowthLawLinear::GrowthLawLinear()
  : GrowthLawDyn(NULL)
{
}


/*----------------------------------------------------------------------------*/
MAT::GrowthLawLinear::GrowthLawLinear(MAT::PAR::GrowthLawLinear* params)
  : GrowthLawDyn(params)
{
}



/*----------------------------------------------------------------------------*/
void MAT::GrowthLawLinear::EvaluateGrowthFunction
(
  double & growthfunc,
  double traceM,
  double theta
)
{
  MAT::PAR::GrowthLawLinear* params = Parameter();
  // parameters
  const double hommandel   = params->hommandel_;
  const double thetaplus   = params->thetaplus_; //1.5;
  const double kthetaplus  = params->kthetaplus_; //1.0; 0.5;
  const double mthetaplus  = params->mthetaplus_; //2.0; 4.0;
  const double thetaminus  = params->thetaminus_; //0.5;
  const double kthetaminus = params->kthetaminus_; //2.0; 0.25;
  const double mthetaminus = params->mthetaminus_; //3.0; 5.0;

  double ktheta  = 0.0;

  if (traceM > hommandel) {
    ktheta  = kthetaplus*pow((thetaplus-theta)/(thetaplus-1.0),mthetaplus);
  } else if (traceM < hommandel) {
    ktheta  = kthetaminus*pow((theta-thetaminus)/(1.0-thetaminus),mthetaminus);
  }

  growthfunc = GetCfac()[0] * ktheta * (traceM - hommandel);
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawLinear::EvaluateGrowthFunctionDerivTheta
(
  double & dgrowthfunctheta,
  double traceM,
  double theta,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
  const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
)
{
  MAT::PAR::GrowthLawLinear* params = Parameter();
  // parameters
  const double hommandel   = params->hommandel_;
  const double thetaplus   = params->thetaplus_;
  const double kthetaplus  = params->kthetaplus_;
  const double mthetaplus  = params->mthetaplus_;
  const double thetaminus  = params->thetaminus_;
  const double kthetaminus = params->kthetaminus_;
  const double mthetaminus = params->mthetaminus_;

  double ktheta  = 0.0;
  double dktheta = 0.0;

  if (traceM > hommandel) {
    ktheta  = kthetaplus*pow((thetaplus-theta)/(thetaplus-1.0),mthetaplus);
    dktheta = mthetaplus*kthetaplus*pow((thetaplus-theta)/(thetaplus-1.0),mthetaplus-1.0)/(1.0-thetaplus);
  } else if (traceM < hommandel) {
    ktheta  = kthetaminus*pow((theta-thetaminus)/(1.0-thetaminus),mthetaminus);
    dktheta = mthetaminus*kthetaminus*pow((theta-thetaminus)/(1.0-thetaminus),mthetaminus-1.0)/(1.0-thetaminus);
  }

  double temp = 0.0;
  for (int i = 0; i < 6; i++)
  {
    temp += Cdach(i) *
            (   cmatelastic(i, 0) * Cdach(0) + cmatelastic(i, 1) * Cdach(1)
              + cmatelastic(i, 2) * Cdach(2) + cmatelastic(i, 3) * Cdach(3)
              + cmatelastic(i, 4) * Cdach(4) + cmatelastic(i, 5) * Cdach(5));
  }

  double dtraceM =  -(2.0 * traceM + temp) / theta;

  dgrowthfunctheta = GetCfac()[0] * (dktheta * (traceM - hommandel) + ktheta * dtraceM);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawLinear::EvaluateGrowthFunctionDerivC
(
  LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
  double traceM,
  double theta,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
  const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
)
{
  MAT::PAR::GrowthLawLinear* params = Parameter();
  // parameters
  const double hommandel   = params->hommandel_;
  const double thetaplus   = params->thetaplus_; //1.5;
  const double kthetaplus  = params->kthetaplus_; //1.0; 0.5;
  const double mthetaplus  = params->mthetaplus_; //2.0; 4.0;
  const double thetaminus  = params->thetaminus_; //0.5;
  const double kthetaminus = params->kthetaminus_; //2.0; 0.25;
  const double mthetaminus = params->mthetaminus_; //3.0; 5.0;

  double ktheta  = 0.0;

  if (traceM > hommandel) {
    ktheta  = kthetaplus*pow((thetaplus-theta)/(thetaplus-1.0),mthetaplus);
  } else if (traceM < hommandel) {
    ktheta  = kthetaminus*pow((theta-thetaminus)/(1.0-thetaminus),mthetaminus);
  }

  for (int j = 0; j < NUM_STRESS_3D; j++)
  {
    double Ccmatelasj =    cmat(0, j) * C(0) + cmat(1, j) * C(1)
                         + cmat(2, j) * C(2) + cmat(3, j) * C(3)
                         + cmat(4, j) * C(4) + cmat(5, j) * C(5);

    dgrowthfuncdC(j) = GetCfac()[0] * ktheta * (S(j)+0.5*Ccmatelasj);
  }

  return;
}

Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawExp::CreateMaterial()
{
  return Teuchos::null;
}

Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawExp::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawExp(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawExp::GrowthLawExp(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: GrowthLawDyn(matdata),
  mandel_(matdata->GetDouble("MANDEL"))
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawExp::GrowthLawExp()
  : GrowthLawDyn(NULL)
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawExp::GrowthLawExp(MAT::PAR::GrowthLawExp* params)
  : GrowthLawDyn(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawExp::EvaluateGrowthFunction
(
  double & growthfunc,
  double traceM,
  double theta
)
{
  MAT::PAR::GrowthLawExp* params = Parameter();
  // parameters
  const double mandel   = params->mandel_;

  growthfunc = GetCfac()[0] * theta * exp(-traceM*traceM/(mandel*mandel));
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawExp::EvaluateGrowthFunctionDerivTheta
(
  double & dgrowthfunctheta,
  double traceM,
  double theta,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
  const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
)
{
  MAT::PAR::GrowthLawExp* params = Parameter();
  // parameters
  const double mandel   = params->mandel_;

  double temp = 0.0;
  for (int i = 0; i < 6; i++)
  {
    temp += Cdach(i) *
            (   cmatelastic(i, 0) * Cdach(0) + cmatelastic(i, 1) * Cdach(1)
              + cmatelastic(i, 2) * Cdach(2) + cmatelastic(i, 3) * Cdach(3)
              + cmatelastic(i, 4) * Cdach(4) + cmatelastic(i, 5) * Cdach(5));
  }

  double dtraceM =  -(2.0 * traceM + temp) / theta;

  dgrowthfunctheta =  GetCfac()[0] * exp(-traceM*traceM/(mandel*mandel)) *
                      (
                            1.0
                          - 2.0 * traceM/(mandel*mandel) * theta * dtraceM
                      );

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawExp::EvaluateGrowthFunctionDerivC
(
  LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
  double traceM,
  double theta,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
  const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
)
{
  MAT::PAR::GrowthLawExp* params = Parameter();
  // parameters
  const double mandel   = params->mandel_;

  {
    double expdC = - 2.0 * theta * traceM/(mandel*mandel)
                                   * exp(-traceM*traceM/(mandel*mandel));
    for (int j = 0; j < NUM_STRESS_3D; j++)
    {
      double Ccmatelasj =    cmat(0, j) * C(0) + cmat(1, j) * C(1)
                          + cmat(2, j) * C(2) + cmat(3, j) * C(3)
                          + cmat(4, j) * C(4) + cmat(5, j) * C(5);

      dgrowthfuncdC(j) = GetCfac()[0] * expdC * (S(j)+0.5*Ccmatelasj);
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawBiofilm::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawBiofilm::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawBiofilm(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawBiofilm::GrowthLawBiofilm(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: GrowthLawDyn(matdata),
  mandel_(matdata->GetDouble("MANDEL"))
{

}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawBiofilm::GrowthLawBiofilm()
  : GrowthLawDyn(NULL)
{
}


/*----------------------------------------------------------------------------*/
MAT::GrowthLawBiofilm::GrowthLawBiofilm(MAT::PAR::GrowthLawBiofilm* params)
  : GrowthLawDyn(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawBiofilm::EvaluateGrowthFunction
(
  double & growthfunc,
  double traceM,
  double theta
)
{
  MAT::PAR::GrowthLawBiofilm* params = Parameter();
  // parameters
  const double mandel   = params->mandel_;

  growthfunc = GetCfac()[0] - abs(theta * exp(-traceM*traceM/(mandel*mandel)));
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawBiofilm::EvaluateGrowthFunctionDerivTheta
(
  double & dgrowthfunctheta,
  double traceM,
  double theta,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
  const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
)
{
  MAT::PAR::GrowthLawBiofilm* params = Parameter();
  // parameters
  const double mandel   = params->mandel_;

  double temp = 0.0;
  for (int i = 0; i < 6; i++)
  {
    temp += Cdach(i) *
            (   cmatelastic(i, 0) * Cdach(0) + cmatelastic(i, 1) * Cdach(1)
              + cmatelastic(i, 2) * Cdach(2) + cmatelastic(i, 3) * Cdach(3)
              + cmatelastic(i, 4) * Cdach(4) + cmatelastic(i, 5) * Cdach(5));
  }

  double dtraceM =  -(2.0 * traceM + temp) / theta;

  dgrowthfunctheta = exp(-traceM*traceM/(mandel*mandel)) *
                      (
                            1.0
                          - 2.0 * traceM/(mandel*mandel) * theta * dtraceM
                      );

  double growthfunc = theta * exp(-traceM*traceM/(mandel*mandel));

  if (growthfunc<1.0e-14) dgrowthfunctheta = 0.0;
  else dgrowthfunctheta = -abs(growthfunc)/growthfunc * dgrowthfunctheta;

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawBiofilm::EvaluateGrowthFunctionDerivC
(
  LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
  double traceM,
  double theta,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
  const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
)
{
  MAT::PAR::GrowthLawBiofilm* params = Parameter();
  // parameters
  const double mandel   = params->mandel_;

  {
    double expdC = - 2.0 * theta * traceM/(mandel*mandel)
                                   * exp(-traceM*traceM/(mandel*mandel));
    for (int j = 0; j < NUM_STRESS_3D; j++)
    {
      double Ccmatelasj =    cmat(0, j) * C(0) + cmat(1, j) * C(1)
                          + cmat(2, j) * C(2) + cmat(3, j) * C(3)
                          + cmat(4, j) * C(4) + cmat(5, j) * C(5);

      dgrowthfuncdC(j) = expdC * (S(j)+0.5*Ccmatelasj);
    }
  }

  double growthfunc = theta * exp(-traceM*traceM/(mandel*mandel));

  if (growthfunc<1.0e-14) dgrowthfuncdC.Scale(0.0);
  else dgrowthfuncdC.Scale(-abs(growthfunc)/growthfunc);

  return;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawConst::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawConst::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawConst(this));
}

/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawConst::GrowthLawConst(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  thetarate_(matdata->GetDouble("THETARATE"))
{

}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawConst::GrowthLawConst()
  : GrowthLawStatic(NULL)
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawConst::GrowthLawConst(MAT::PAR::GrowthLawConst* params)
  : GrowthLawStatic(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawConst::Evaluate(double* theta,
                                const double& thetaold,
                                LINALG::Matrix<6,1>* dthetadC,
                                Teuchos::RCP<MAT::So3Material> matelastic,
                                const LINALG::Matrix<3,3>* defgrd,
                                const LINALG::Matrix<6,1>* glstrain,
                                Teuchos::ParameterList& params,
                                const int eleGID)
{
  double dt = params.get<double>("delta time", -1.0);
  *theta = thetaold + Parameter()->thetarate_*dt;

  dthetadC->Scale(0.0);

  return;
}


/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawAC::GrowthLawAC(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  Sc1_(matdata->GetInt("SCALAR1")),
  alpha_(matdata->GetDouble("ALPHA")),
  Sc2_(matdata->GetInt("SCALAR2")),
  beta_(matdata->GetDouble("BETA"))
{
  if (Sc1_ < 1)
    dserror("At least on scalar field must induce growth");
  if (alpha_ < 0)
    dserror("The influence of scalar field SCALAR1 to growth can't be negativ");
  if (Sc2_ < 1)
    dserror("If you choose a second scalar field to induce growth, choose a existing one!!");
  if (beta_ < 0)
    dserror("The influence of scalar field SCALAR2 to growth can't be negativ");
}

Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawAC::CreateMaterial()
{
  return Teuchos::null;
}

Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawAC::CreateGrowthLaw()
{

  return Teuchos::rcp(new MAT::GrowthLawAC(this));
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawAC::GrowthLawAC()
  : GrowthLawStatic(NULL)
{
}


/*----------------------------------------------------------------------------*/
MAT::GrowthLawAC::GrowthLawAC(MAT::PAR::GrowthLawAC* params)
  : GrowthLawStatic(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawAC::Evaluate(double* theta,
                                const double& thetaold,
                                LINALG::Matrix<6,1>* dthetadC,
                                Teuchos::RCP<MAT::So3Material> matelastic,
                                const LINALG::Matrix<3,3>* defgrd,
                                const LINALG::Matrix<6,1>* glstrain,
                                Teuchos::ParameterList& params,
                                const int eleGID)
{
  int Sc1 = Parameter()->Sc1_;
  double alpha= Parameter()->alpha_;
  int Sc2 = Parameter()->Sc2_; //if no second scalar is chosen to induce growth, this
                               //points out the first scalar field
  double beta= Parameter()->beta_; //if no second scalar is chosen to induce growth,
                                   //beta is zero and hence has no influence on growth

  double J = defgrd->Determinant();
  *theta = pow(1.0 + alpha * GetCfac().at(Sc1-1) * J + beta * GetCfac().at(Sc2-1) * J ,0.33333333333333333);

  double tmp1 = (alpha*GetCfac().at(Sc1-1) + beta * GetCfac().at(Sc2-1))*J;
  double tmp2 = tmp1/6.0*pow(1+tmp1,-0.66666666666666666);

  LINALG::Matrix<3, 3> C(true);
  C.MultiplyTN(*defgrd,*defgrd);

  // linearization of growth law
  LINALG::Matrix<3, 3> dThetadC(true);
  dThetadC.Invert(C);
  dThetadC.Scale(tmp2);

  //transform dThetadC into a vector
  (*dthetadC)(0)=dThetadC(0,0);
  (*dthetadC)(1)=dThetadC(1,1);
  (*dthetadC)(2)=dThetadC(2,2);
  (*dthetadC)(3)=2*dThetadC(0,1);
  (*dthetadC)(4)=2*dThetadC(1,2);
  (*dthetadC)(5)=2*dThetadC(0,2);
}


/*----------------------------------------------------------------------------*/
MAT::PAR::GrowthLawACRadial::GrowthLawACRadial(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
  : GrowthLawAC(matdata)
{
}

Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawACRadial::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawACRadial(this));
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawACRadial::GrowthLawACRadial()
  : GrowthLawStatic(NULL)
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthLawACRadial::GrowthLawACRadial(MAT::PAR::GrowthLawAC* params)
  : GrowthLawStatic(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthLawACRadial::Evaluate(double* theta,
                                      const double& thetaold,
                                      LINALG::Matrix<6,1>* dthetadC,
                                      Teuchos::RCP<MAT::So3Material> matelastic,
                                      const LINALG::Matrix<3,3>* defgrd,
                                      const LINALG::Matrix<6,1>* glstrain,
                                      Teuchos::ParameterList& params,
                                      const int eleGID)
{
  int Sc1 = Parameter()->Sc1_;
  double alpha= Parameter()->alpha_;
  int Sc2 = Parameter()->Sc2_; //if no second scalar is chosen to induce growth, this
                               //points out the first scalar field
  double beta= Parameter()->beta_; //if no second scalar is chosen to induce growth,
                                   //beta is zero and hence has no influence on growth

  double J = 1.0; //defgrd->Determinant();

  *theta = 1.0 + alpha * GetCfac().at(Sc1-1) * J + beta * GetCfac().at(Sc2-1) * J;

  double tmp1 = 0.5 * (alpha*GetCfac().at(Sc1-1) + beta * GetCfac().at(Sc2-1))*J;


  LINALG::Matrix<3, 3> C(true);
  C.MultiplyTN(*defgrd,*defgrd);

  // linearization of growth law
  LINALG::Matrix<3, 3> dThetadC(true);
  dThetadC.Invert(C);
  dThetadC.Scale(tmp1);

  //transform dThetadC into a vector
  (*dthetadC)(0)=dThetadC(0,0);
  (*dthetadC)(1)=dThetadC(1,1);
  (*dthetadC)(2)=dThetadC(2,2);
  (*dthetadC)(3)=2*dThetadC(0,1);
  (*dthetadC)(4)=2*dThetadC(1,2);
  (*dthetadC)(5)=2*dThetadC(0,2);

}
