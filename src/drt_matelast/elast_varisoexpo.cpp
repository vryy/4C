/*----------------------------------------------------------------------*/
/*!
\file elast_varisoexpo.cpp
\brief


the input line should read
  MAT 1 ELAST_VarIsoExpo FRAC 0.5 C 100

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_varisoexpo.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::VarIsoExpo::VarIsoExpo(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  frac_(matdata->GetDouble("FRAC")),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::VarIsoExpo::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::VarIsoExpo( this ) );
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VarIsoExpo::VarIsoExpo()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VarIsoExpo::VarIsoExpo(MAT::ELASTIC::PAR::VarIsoExpo* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VarIsoExpo::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  havecoefficients = havecoefficients or true;

  const double frac = params_ -> frac_;
  const double k1 = params_ -> k1_;
  const double k2 = params_ -> k2_;

  gamma(0) += (2.*modinv(0)*frac*k1-6.*k2)*exp(k2*(9.-6*modinv(0)+pow(modinv(0), 2)));
  delta(0) += 4.*frac*k1*(1.+18.*k2+2.*k2*pow(modinv(0), 2)-12.*k2*modinv(0))* exp(k2*(9.-6*modinv(0)+pow(modinv(0), 2)));

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
