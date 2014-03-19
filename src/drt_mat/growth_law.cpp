/*----------------------------------------------------------------------*/
/*!
 \file growth_law.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#include "growth_law.H"
#include "matpar_material.H"

Teuchos::RCP<MAT::Material> MAT::PAR::GrowthLawLinear::CreateMaterial()
{
  return Teuchos::null;
}

Teuchos::RCP<MAT::GrowthLaw> MAT::PAR::GrowthLawLinear::CreateGrowthLaw()
{
  return Teuchos::rcp(new MAT::GrowthLawLinear(this));
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::GrowthLawLinear::GrowthLawLinear(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  thetaplus_(matdata->GetDouble("THETAPLUS")),
  kthetaplus_(matdata->GetDouble("KPLUS")),
  mthetaplus_(matdata->GetDouble("MPLUS")),
  thetaminus_(matdata->GetDouble("THETAMINUS")),
  kthetaminus_(matdata->GetDouble("KMINUS")),
  mthetaminus_(matdata->GetDouble("MMINUS")),
  hommandel_(matdata->GetDouble("HOMMANDEL"))
{

}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawLinear::GrowthLawLinear()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawLinear::GrowthLawLinear(MAT::PAR::GrowthLawLinear* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 |  Evaluate growth function                           (protected)        02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthLawLinear::EvaluateGrowthFunction
(
  double & growthfunc,
  double traceM,
  double theta
)
{
  // parameters
  const double hommandel   = params_->hommandel_;
  const double thetaplus   = params_->thetaplus_; //1.5;
  const double kthetaplus  = params_->kthetaplus_; //1.0; 0.5;
  const double mthetaplus  = params_->mthetaplus_; //2.0; 4.0;
  const double thetaminus  = params_->thetaminus_; //0.5;
  const double kthetaminus = params_->kthetaminus_; //2.0; 0.25;
  const double mthetaminus = params_->mthetaminus_; //3.0; 5.0;

  double ktheta  = 0.0;

  if (traceM > hommandel) {
    ktheta  = kthetaplus*pow((thetaplus-theta)/(thetaplus-1.0),mthetaplus);
  } else if (traceM < hommandel) {
    ktheta  = kthetaminus*pow((theta-thetaminus)/(1.0-thetaminus),mthetaminus);
  }

  growthfunc = ktheta * (traceM - hommandel);
}

/*----------------------------------------------------------------------*
 |  Evaluate derivative of growth function            (protected)   02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthLawLinear::EvaluateGrowthFunctionDerivTheta
(
  double & dgrowthfunctheta,
  double traceM,
  double theta,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
  const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
)
{
  // parameters
  const double hommandel   = params_->hommandel_;
  const double thetaplus   = params_->thetaplus_; //1.5;
  const double kthetaplus  = params_->kthetaplus_; //1.0; 0.5;
  const double mthetaplus  = params_->mthetaplus_; //2.0; 4.0;
  const double thetaminus  = params_->thetaminus_; //0.5;
  const double kthetaminus = params_->kthetaminus_; //2.0; 0.25;
  const double mthetaminus = params_->mthetaminus_; //3.0; 5.0;

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

  dgrowthfunctheta = dktheta * (traceM - hommandel) + ktheta * dtraceM;

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate derivative of growth function            (protected)   02/10|
 *----------------------------------------------------------------------*/
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
  // parameters
  const double hommandel   = params_->hommandel_;
  const double thetaplus   = params_->thetaplus_; //1.5;
  const double kthetaplus  = params_->kthetaplus_; //1.0; 0.5;
  const double mthetaplus  = params_->mthetaplus_; //2.0; 4.0;
  const double thetaminus  = params_->thetaminus_; //0.5;
  const double kthetaminus = params_->kthetaminus_; //2.0; 0.25;
  const double mthetaminus = params_->mthetaminus_; //3.0; 5.0;

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

    dgrowthfuncdC(j) = ktheta * (S(j)+0.5*Ccmatelasj);
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


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::GrowthLawExp::GrowthLawExp(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  mandel_(matdata->GetDouble("MANDEL"))
{

}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawExp::GrowthLawExp()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawExp::GrowthLawExp(MAT::PAR::GrowthLawExp* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 |  Evaluate growth function                           (protected)        02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthLawExp::EvaluateGrowthFunction
(
  double & growthfunc,
  double traceM,
  double theta
)
{
  // parameters
  const double mandel   = params_->mandel_;

  growthfunc = theta * exp(-traceM*traceM/(mandel*mandel));
}

/*----------------------------------------------------------------------*
 |  Evaluate derivative of growth function            (protected)   02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthLawExp::EvaluateGrowthFunctionDerivTheta
(
  double & dgrowthfunctheta,
  double traceM,
  double theta,
  const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
  const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
)
{
  // parameters
  const double mandel   = params_->mandel_;

  double temp = 0.0;
  for (int i = 0; i < 6; i++)
  {
    temp += Cdach(i) *
            (   cmatelastic(i, 0) * Cdach(0) + cmatelastic(i, 1) * Cdach(1)
              + cmatelastic(i, 2) * Cdach(2) + cmatelastic(i, 3) * Cdach(3)
              + cmatelastic(i, 4) * Cdach(4) + cmatelastic(i, 5) * Cdach(5));
  }

  double dtraceM =  -(2.0 * traceM + temp) / theta;

  dgrowthfunctheta =  exp(-traceM*traceM/(mandel*mandel)) *
                      (
                            1.0
                          - 2.0 * traceM/(mandel*mandel) * theta * dtraceM
                      );

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate derivative of growth function            (protected)   02/10|
 *----------------------------------------------------------------------*/
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
  // parameters
  const double mandel   = params_->mandel_;

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

  return;
}
