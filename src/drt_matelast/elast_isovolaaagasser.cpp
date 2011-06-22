/*----------------------------------------------------------------------*/
/*!
\file elast_isovolaaagasser.cpp
\brief

This file contains the routines required to calculate the isochoric contribution
of the aaagasser material and the corresponding volumetric contribution.

MAT 20 ELAST_IsoVolAAAGasser CLUM 2.62E3 CMED 1.98E3 CABLUM 1.73E3 NUE 0.49 BETA -2.0

<pre>
maintainer: Andreas Maier
            a.maier@lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isovolaaagasser.H"


/*----------------------------------------------------------------------*
 |                                                         AMaier  06/11|
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoVolAAAGasser::IsoVolAAAGasser(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  clum_(matdata->GetDouble("CLUM")),
  cmed_(matdata->GetDouble("CMED")),
  cablum_(matdata->GetDouble("CABLUM")),
  nue_(matdata->GetDouble("NUE")),
  beta_(matdata->GetDouble("BETA"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::IsoVolAAAGasser::CreateMaterial()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public) AMaier 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVolAAAGasser::IsoVolAAAGasser()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                 (public)   AMaier 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVolAAAGasser::IsoVolAAAGasser(MAT::ELASTIC::PAR::IsoVolAAAGasser* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolAAAGasser::AddCoefficientsPrincipal(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv,
  double normdist
  )
{
  havecoefficients = havecoefficients or true;
  
  if (normdist==-999.0) dserror("Aneurysm mean ilt distance not found. Did you switch on 'PATSPEC'? (besides other possible errors of course)");

  double cele;
  if (0.0 <= normdist <= 0.5)
    cele = ((normdist - 0.5)/(-0.5))* params_->clum_ + (normdist/0.5)*params_->cmed_;
  else if(0.5 < normdist <= 1.0)
    cele = ((normdist - 1.0)/(-0.5))*params_->cmed_ + ((normdist - 0.5)/0.5)*params_->cablum_;
  else
    dserror("Unable to calculate valid stiffness parameter in material AAAGasser");

  // principal coefficients
  gamma(0) += 0.;
  gamma(1) += cele * 4.*pow(prinv(2),-2./3.);
  gamma(2) += cele * (-(4./3.)*pow(prinv(2),-2./3.)*(prinv(0)*prinv(0) - 2.*prinv(1)) + (2./(1.-2.*params_->nue_))*4.*(1.-pow(prinv(2),-params_->beta_/2.))/params_->beta_);
  delta(0) += 0.0;
  delta(1) += 0.0;
  delta(2) += 0.0;
  delta(3) += 0.0;
  delta(4) += cele * (-16./3.)*pow(prinv(2),-2./3.);
  delta(5) += cele * ((16./9.)*pow(prinv(2),-2./3.)*(prinv(0)*prinv(0) - 2.*prinv(1)) + (2./(1.-2.*params_->nue_))*4.*pow(prinv(2),-params_->beta_/2.));
  delta(6) += cele * ((8./3.)*pow(prinv(2),-2./3.)*(prinv(0)*prinv(0) - 2.*prinv(1)) + (2./(1.-2.*params_->nue_))*4.*2.*(pow(prinv(2),-params_->beta_/2.)-1.)/params_->beta_);
  delta(7) += cele * 8.*pow(prinv(2),-2./3.);

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
