/*----------------------------------------------------------------------*/
/*!
\file elast_vologden.cpp
\brief


the input line should read
  MAT 1 ELAST_VolOgden KAPPA 100 BETA -2

<pre>
Maintainer: Sophie Rausch
            rausch,kloeppel@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_vologden.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELAST::PAR::VolOgden::VolOgden(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  kappa_(matdata->GetDouble("KAPPA")),
  beta_(matdata->GetDouble("BETA"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELAST::VolOgden::VolOgden()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELAST::VolOgden::VolOgden(MAT::ELAST::PAR::VolOgden* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELAST::VolOgden::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  havecoefficients = havecoefficients or true;

  const double kappa = params_ -> kappa_;
  const double beta = params_ -> beta_;

  gamma(2) += kappa*(1-pow(modinv(2), -beta))/(modinv(2)*beta);
  delta(4) += kappa*pow(modinv(2), -beta-1);

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
