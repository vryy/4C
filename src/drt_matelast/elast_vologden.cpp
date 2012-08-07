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

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_vologden.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::VolOgden::VolOgden(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  kappa_(matdata->GetDouble("KAPPA")),
  beta_(matdata->GetDouble("BETA"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::VolOgden::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::VolOgden( this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VolOgden::VolOgden()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VolOgden::VolOgden(MAT::ELASTIC::PAR::VolOgden* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolOgden::AddCoefficientsModified(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{

  const double kappa = params_ -> kappa_;
  const double beta = params_ -> beta_;

  gamma(2) += kappa*(1-pow(modinv(2), -beta))/(modinv(2)*beta);
  delta(4) += kappa*pow(modinv(2), -beta-1);

  return;
}


/*----------------------------------------------------------------------*/
