/*----------------------------------------------------------------------*/
/*!
\file elast_isoexpo.cpp
\brief

This file contains the routines required to calculate the isochoric contribution
of an exponential  material

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isoexpo.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoExpo::IsoExpo(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoExpo::IsoExpo()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoExpo::IsoExpo(MAT::ELASTIC::PAR::IsoExpo* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoExpo::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  havecoefficients = havecoefficients or true;

  const double k1 = params_ -> k1_;
  const double k2 = params_ -> k2_;

  gamma(0) += (2.*modinv(0)*k1-6.*k2)*exp(k2*(9.-6*modinv(0)+pow(modinv(0), 2)));

  delta(0) += 4.*k1*(1.+18.*k2+2.*k2*pow(modinv(0), 2)-12.*k2*modinv(0))* exp(k2*(9.-6*modinv(0)+pow(modinv(0), 2)));

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
