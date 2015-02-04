/*----------------------------------------------------------------------*/
/*!
 \file growth_law.cpp

 \brief
This file contains routines needed for the calculation of the volumetric growth
parameter theta. The returned results can be either theta itself or its derivative
 <pre>
   Maintainer: Moritz Thon
               thon@lnm.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#include "growth_law.H"
#include "matpar_material.H"


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLaw::GrowthLaw()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLaw::GrowthLaw(MAT::PAR::Parameter* params)
  : params_(params)
{
}

/// This is a law returning the derivatives of theta
/// See e.g.:
///- Lubarda, V. & Hoger, A., On the mechanics of solids with a growing mass,
///  International Journal of Solids and Structures, 2002, 39, 4627-4664
///- Himpel, G.; Kuhl, E.; Menzel, A. & Steinmann, P., Computational modelling
///  of isotropic multiplicative growth, Computer Modeling in Engineering
///  and Sciences, 2005, 8, 119-134

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
  : GrowthLaw(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawLinear::GrowthLawLinear(MAT::PAR::GrowthLawLinear* params)
  : GrowthLaw(params)
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

    dgrowthfuncdC(j) = ktheta * (S(j)+0.5*Ccmatelasj);
  }

  return;
}

/// This is a law returning the derivatives of theta
/// See e.g.: ???
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
  : GrowthLaw(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawExp::GrowthLawExp(MAT::PAR::GrowthLawExp* params)
  : GrowthLaw(params)
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
  MAT::PAR::GrowthLawExp* params = Parameter();
  // parameters
  const double mandel   = params->mandel_;

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

      dgrowthfuncdC(j) = expdC * (S(j)+0.5*Ccmatelasj);
    }
  }

  return;
}

/// This is a growth law directly calculating theta itself

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawAC::GrowthLawAC()
  : GrowthLaw(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawAC::GrowthLawAC(MAT::PAR::GrowthLawAC* params)
  : GrowthLaw(params)
{
}

/*----------------------------------------------------------------------*
 |Calculate Theta from the given mean concentrations          thon 11/14|
 *----------------------------------------------------------------------*/
double MAT::GrowthLawAC::CalculateTheta
(
    Teuchos::RCP<std::vector<double> > concentrations, //pointer to a vector containing element wise mean concentrations
    const double J
)
{
  int Sc1 = Parameter()->Sc1_;
  double alpha= Parameter()->alpha_;
  int Sc2 = Parameter()->Sc2_; //if no second scalar is chosen to induce growth, this points out the first scalar field
  double beta= Parameter()->beta_; //if no second scalar is chosen to induce growth, beta is zero and hence has no influence on growth

  double theta = pow(1.0 + alpha * concentrations->at(Sc1-1) * J + beta * concentrations->at(Sc2-1) * J ,0.33333333333333333);

  return theta;
}

/*----------------------------------------------------------------------*
 | Calculate derivative of Theta w.r.t right cauchy-green strain   thon 01/15|
 *----------------------------------------------------------------------*/
void MAT::GrowthLawAC::CalculateThetaDerivC
(
    LINALG::Matrix<3,3>& dThetadC,
    const LINALG::Matrix<3,3>& C, ///< cauchy-green strains
    Teuchos::RCP<std::vector<double> > concentrations, ///< mean concentrations
    const double J ///< det(F)
)
{
  // dTheta/dC = 1/6 * \alpha*c(i)*J * (1+\alpha*c(i)*J)^(-2/3) * C^(-1)
  int Sc1 = Parameter()->Sc1_;
  double alpha= Parameter()->alpha_;
  int Sc2 = Parameter()->Sc2_; //if no second scalar is chosen to induce growth, this points out the first scalar field
  double beta= Parameter()->beta_; //if no second scalar is chosen to induce growth, beta is zero and hence has no influence on growth

  double tmp1 = (alpha*concentrations->at(Sc1-1) + beta * concentrations->at(Sc2-1))*J;
  double tmp2 = tmp1/6.0*pow(1+tmp1,-0.66666666666666666);

  dThetadC.Invert(C);
  dThetadC.Scale(tmp2);
}


/// This is a growth law directly calculating theta itself

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawACRadial::GrowthLawACRadial()
  : GrowthLaw(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthLawACRadial::GrowthLawACRadial(MAT::PAR::GrowthLawAC* params)
  : GrowthLaw(params)
{
}

/*----------------------------------------------------------------------*
 |Calculate Theta from the given mean concentrations          thon 11/14|
 *----------------------------------------------------------------------*/
double MAT::GrowthLawACRadial::CalculateTheta
(
    Teuchos::RCP<std::vector<double> > concentrations, //pointer to a vector containing element wise mean concentrations
    const double J
)
{
  // Theta = 1 + \alpha*c(i)* J
  int Sc1 = Parameter()->Sc1_;
  double alpha= Parameter()->alpha_;
  int Sc2 = Parameter()->Sc2_; //if no second scalar is chosen to induce growth, this points out the first scalar field
  double beta= Parameter()->beta_; //if no second scalar is chosen to induce growth, beta is zero and hence has no influence on growth

  double theta = 1.0 + alpha * concentrations->at(Sc1-1) * J + beta * concentrations->at(Sc2-1) * J;

  return theta;
}

/*----------------------------------------------------------------------*
 | Calculate derivative of Theta w.r.t right cauchy-green strain   thon 01/15|
 *----------------------------------------------------------------------*/
void MAT::GrowthLawACRadial::CalculateThetaDerivC
(
    LINALG::Matrix<3,3>& dThetadC,
    const LINALG::Matrix<3,3>& C, ///< cauchy-green strain
    Teuchos::RCP<std::vector<double> > concentrations, ///< mean concentrations
    const double J ///< det(F)
)
{
  // dTheta/dC = 1/2 * \alpha*c(i)*J * C^(-1)
  int Sc1 = Parameter()->Sc1_;
  double alpha= Parameter()->alpha_;
  int Sc2 = Parameter()->Sc2_; //if no second scalar is chosen to induce growth, this points out the first scalar field
  double beta= Parameter()->beta_; //if no second scalar is chosen to induce growth, beta is zero and hence has no influence on growth

  double tmp1 = 0.5 * (alpha*concentrations->at(Sc1-1) + beta * concentrations->at(Sc2-1))*J;

  dThetadC.Invert(C);
  dThetadC.Scale(tmp1);
}
