/*----------------------------------------------------------------------*/
/*!
\file elast_isoexpopow.cpp
\brief

This file contains the routines required to calculate the isochoric contribution
of an exponential  material

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isoexpopow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoExpoPow::IsoExpoPow(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2")),
  d_(matdata->GetInt("C"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::IsoExpoPow::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::IsoExpoPow( this ) );
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoExpoPow::IsoExpoPow()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoExpoPow::IsoExpoPow(MAT::ELASTIC::PAR::IsoExpoPow* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoExpoPow::AddCoefficientsModified(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{

  const double k1 = params_ -> k1_;
  const double k2 = params_ -> k2_;
  const int d = params_->d_;
  if (d == 0) dserror("D=0 is not valid, choose a bigger D");

  //gamma(0) += (2.*modinv(0)*k1-6.*k1)*exp(k2*(9.-6*modinv(0)+pow(modinv(0), 2)));
  gamma(0) += k1*d* pow(modinv(0)-3.0,d-1)* exp(k2*pow(modinv(0)-3.,d));

  //delta(0) += 4.*k1*(1.+18.*k2+2.*k2*pow(modinv(0), 2)-12.*k2*modinv(0))* exp(k2*(9.-6*modinv(0)+pow(modinv(0), 2)));
  if (d >= 2)
    delta(0) += 2.*k1*d*((d-1)*pow(modinv(0)-3.,d-2) + k2*d*pow(modinv(0)-3.0,2*d-2))* exp(k2*pow(modinv(0)-3.,d));
  else
    delta(0) += 2.*k1*d*d*k2* exp(k2*pow(modinv(0)-3.,d));

  return;
}
