/*----------------------------------------------------------------------*/
/*!
 \file poro_law.cpp

 \brief calculation classes for evaluation of constitutive relation for porosity

 <pre>
 \level 2

 \maintainer Anh-Tu Vuong
             vuong@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#include "matpar_bundle.H"
#include "poro_law.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroLaw::PoroLaw(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroLawLinear::PoroLawLinear(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: PoroLaw(matdata),
  bulkmodulus_(matdata->GetDouble("BULKMODULUS"))
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PoroLawLinear::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawLinear::ComputePorosity(
                                       const double& refporosity,
                                       const double& press,
                                       const double& J,
                                       const int& gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       double* dphi_dphiref)
{
  porosity = refporosity + 1.0 / bulkmodulus_ * press
                + (1.0-refporosity) * (J-1.0);

  if(dphi_dp)      *dphi_dp      = 1.0 / bulkmodulus_;
  if(dphi_dJ)      *dphi_dJ      = (1.0-refporosity);
  if(dphi_dJdp)    *dphi_dJdp    = 0.0;
  if(dphi_dJJ)     *dphi_dJJ     = 0.0;
  if(dphi_dpp)     *dphi_dpp     = 0.0;
  if(dphi_dphiref) *dphi_dphiref = 2.0-J;

  return;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawLinear::ConstitutiveDerivatives(
    Teuchos::ParameterList& params,
    double     press,
    double     J,
    double     porosity,
    double     refporosity,
    double*    dW_dp,
    double*    dW_dphi,
    double*    dW_dJ,
    double*    dW_dphiref,
    double*    W)
{
  if(W)            *W            = bulkmodulus_ * (porosity - refporosity - (1.0-refporosity) * (J-1.0))-press;
  if(dW_dp)        *dW_dp        = -1.0;
  if(dW_dphi)      *dW_dphi      = bulkmodulus_;
  if(dW_dJ)        *dW_dJ        = -bulkmodulus_*(1.0-refporosity);
  if(dW_dphiref)   *dW_dphiref   = bulkmodulus_*(-2.0+J);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroLawNeoHooke::PoroLawNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: PoroLaw(matdata),
  bulkmodulus_(matdata->GetDouble("BULKMODULUS")),
  penaltyparameter_(matdata->GetDouble("PENALTYPARAMETER"))
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PoroLawNeoHooke::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawNeoHooke::ComputePorosity(
                                       const double& refporosity,
                                       const double& press,
                                       const double& J,
                                       const int& gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       double* dphi_dphiref)
{
  const double & bulkmodulus  = bulkmodulus_;
  const double & penalty      = penaltyparameter_;

  const double a = (bulkmodulus / (1 - refporosity) + press - penalty / refporosity) * J;
  const double b = -a + bulkmodulus + penalty;
  const double c = b * b  + 4.0 * penalty * a;
  double d = sqrt(c);


  double test = 1 / (2.0 * a) * (-b + d);
  double sign = 1.0;
  if (test >= 1.0 or test < 0.0)
  {
    sign = -1.0;
    d = sign * d;
  }

  const double a_inv = 1.0/a;
  const double d_inv = 1.0/d;
  const double J_inv = 1.0/J;

  const double phi = 1 / (2 * a) * (-b + d);


  const double d_p = J * (-b+2.0*penalty) * d_inv;
  const double d_p_p = ( d * J + d_p * (b - 2.0*penalty) ) * d_inv * d_inv * J;
  const double d_J = a * J_inv * ( -b + 2.0*penalty ) * d_inv;
  const double d_J_p = (d_p * J_inv + ( 1-d_p*d_p*J_inv*J_inv ) *d_inv *a);
  const double d_J_J = ( a*a*J_inv*J_inv-d_J*d_J )* d_inv;

  //d(porosity) / d(p)
  if(dphi_dp) *dphi_dp = (- J * phi + 0.5*(J+d_p))*a_inv;

  //d(porosity) / d(J)
  if(dphi_dJ) *dphi_dJ= (-phi+ 0.5) * J_inv + 0.5*d_J * a_inv;

  //d(porosity) / d(J)d(pressure)
  if(dphi_dJdp) *dphi_dJdp= -J_inv* (*dphi_dp)+ 0.5 * d_J_p * a_inv - 0.5 * d_J*J* a_inv* a_inv ;

  //d^2(porosity) / d(J)^2
  if(dphi_dJJ) *dphi_dJJ= phi*J_inv*J_inv - (*dphi_dJ)*J_inv - 0.5*J_inv*J_inv - 0.5*d_J*J_inv*a_inv + 0.5*d_J_J*a_inv;

  //d^2(porosity) / d(pressure)^2
  if(dphi_dpp) *dphi_dpp= -J*a_inv* (*dphi_dp) + phi*J*J*a_inv*a_inv - 0.5*J*a_inv*a_inv*(J+d_p) + 0.5*d_p_p*a_inv;

  porosity= phi;

  if(dphi_dphiref)
  {
    const double dadphiref = J*(bulkmodulus / ((1 - refporosity)*(1 - refporosity)) + penalty / (refporosity*refporosity));
    const double tmp = 2*dadphiref*a_inv * (-b*(a+b)*a_inv - 2*penalty);
    const double dddphiref = sign*(dadphiref * sqrt(c)*a_inv + tmp);

    *dphi_dphiref = ( a * (dadphiref+dddphiref) - dadphiref * (-b + d) )*0.5*a_inv*a_inv;
  }

  return;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawNeoHooke::ConstitutiveDerivatives(
    Teuchos::ParameterList& params,
    double     press,
    double     J,
    double     porosity,
    double     refporosity,
    double*    dW_dp,
    double*    dW_dphi,
    double*    dW_dJ,
    double*    dW_dphiref,
    double*    W)
{
  //some intermediate values
  const double a = bulkmodulus_ / (1 - refporosity) + press - penaltyparameter_ / refporosity;
  const double b = -1.0*J*a+bulkmodulus_+penaltyparameter_;

  const double scale = 1.0/bulkmodulus_;

  //scale everything with 1/bulkmodulus (I hope this will help the solver...)
  if(W)       *W       = (J*a*porosity*porosity + porosity* b - penaltyparameter_) * scale;
  if(dW_dp)   *dW_dp   = (-1.0*J*porosity *(1.0-porosity)) * scale;
  if(dW_dphi) *dW_dphi = (2.0*J*a*porosity + b) * scale;
  if(dW_dJ)   *dW_dJ   = (a*porosity*porosity - porosity*a) * scale;

  if(dW_dphiref)
  {
    const double dadphiref = J*(bulkmodulus_ / ((1 - refporosity)*(1 - refporosity)) + penaltyparameter_ / (refporosity*refporosity));
    const double dbdphiref = -1.0*J*dadphiref;

    *dW_dphiref = (J*dadphiref*porosity*porosity + porosity* dbdphiref ) * scale;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroLawConstant::PoroLawConstant(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: PoroLaw(matdata)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PoroLawConstant::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawConstant::ComputePorosity(
                                       const double& refporosity,
                                       const double& press,
                                       const double& J,
                                       const int& gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       double* dphi_dphiref)
{
  //porosity is constant -> derivates are zero
  porosity = refporosity;
  if(dphi_dp)      *dphi_dp      = 0.0;
  if(dphi_dJ)      *dphi_dJ      = 0.0;
  if(dphi_dJdp)    *dphi_dJdp    = 0.0;
  if(dphi_dJJ)     *dphi_dJJ     = 0.0;
  if(dphi_dpp)     *dphi_dpp     = 0.0;
  if(dphi_dphiref) *dphi_dphiref = 1.0;

  return;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawConstant::ConstitutiveDerivatives(
    Teuchos::ParameterList& params,
    double     press,
    double     J,
    double     porosity,
    double     refporosity,
    double*    dW_dp,
    double*    dW_dphi,
    double*    dW_dJ,
    double*    dW_dphiref,
    double*    W)
{
  if(W)           *W            = porosity-refporosity;
  if(dW_dp)       *dW_dp        = 0.0;
  if(dW_dphi)     *dW_dphi      = 1.0;
  if(dW_dJ)       *dW_dJ        = 0.0;
  if(dW_dphiref)  *dW_dphiref   = -1.0;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroLawIncompSkeleton::PoroLawIncompSkeleton(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: PoroLaw(matdata)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PoroLawIncompSkeleton::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawIncompSkeleton::ComputePorosity(
                                       const double& refporosity,
                                       const double& press,
                                       const double& J,
                                       const int& gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       double* dphi_dphiref)
{
  porosity = 1.0 - (1.0-refporosity)/J;

  if(dphi_dp)      *dphi_dp      = 0.0;
  if(dphi_dJ)      *dphi_dJ      = (1.0-refporosity)/(J*J);
  if(dphi_dJdp)    *dphi_dJdp    = 0.0;
  if(dphi_dJJ)     *dphi_dJJ     = -2.0*(1.0-refporosity)/(J*J*J);
  if(dphi_dpp)     *dphi_dpp     = 0.0;
  if(dphi_dphiref) *dphi_dphiref = 1.0/J;

  return;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawIncompSkeleton::ConstitutiveDerivatives(
    Teuchos::ParameterList& params,
    double     press,
    double     J,
    double     porosity,
    double     refporosity,
    double*    dW_dp,
    double*    dW_dphi,
    double*    dW_dJ,
    double*    dW_dphiref,
    double*    W)
{
  if(W)           *W            = J*(1.0-porosity)-(1.0-refporosity);
  if(dW_dp)       *dW_dp        = 0.0;
  if(dW_dphi)     *dW_dphi      = -1.0*J;
  if(dW_dJ)       *dW_dJ        = 1.0-porosity;
  if(dW_dphiref)  *dW_dphiref   = 1.0;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroLawLinBiot::PoroLawLinBiot(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: PoroLaw(matdata),
  invBiotModulus_(matdata->GetDouble("INVBIOTMODULUS")),
  biotCoeff_(matdata->GetDouble("BIOTCEOFF"))
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PoroLawLinBiot::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawLinBiot::ComputePorosity(
                                       const double& refporosity,
                                       const double& press,
                                       const double& J,
                                       const int& gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       double* dphi_dphiref)
{
  porosity = refporosity + invBiotModulus_*press + biotCoeff_*(J-1);

  if(dphi_dp)      *dphi_dp      = invBiotModulus_;
  if(dphi_dJ)      *dphi_dJ      = biotCoeff_;
  if(dphi_dJdp)    *dphi_dJdp    = 0.0;
  if(dphi_dJJ)     *dphi_dJJ     = 0.0;
  if(dphi_dpp)     *dphi_dpp     = 0.0;
  if(dphi_dphiref) *dphi_dphiref = 1.0;

  return;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void MAT::PAR::PoroLawLinBiot::ConstitutiveDerivatives(
    Teuchos::ParameterList& params,
    double     press,
    double     J,
    double     porosity,
    double     refporosity,
    double*    dW_dp,
    double*    dW_dphi,
    double*    dW_dJ,
    double*    dW_dphiref,
    double*    W)
{
  if(W)           *W            = porosity - refporosity - invBiotModulus_*press - biotCoeff_*(J-1);
  if(dW_dp)       *dW_dp        = -1.0*invBiotModulus_;
  if(dW_dphi)     *dW_dphi      = 1.0;
  if(dW_dJ)       *dW_dJ        = -1.0*biotCoeff_;
  if(dW_dphiref)  *dW_dphiref   = -1.0;

  return;
}
