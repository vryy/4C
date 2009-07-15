/*----------------------------------------------------------------------*/
/*!
\file elast_volsussmanbathe.cpp
\brief


the input line should read
  MAT 1 ELAST_VolSussmanBathe KAPPA 100

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            rausch,kloeppel@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_volsussmanbathe.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::VolSussmanBathe::VolSussmanBathe(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  kappa_(matdata->GetDouble("KAPPA"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VolSussmanBathe::VolSussmanBathe()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VolSussmanBathe::VolSussmanBathe(MAT::ELASTIC::PAR::VolSussmanBathe* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolSussmanBathe::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  havecoefficients = havecoefficients or true;

  const double kappa = params_ -> kappa_;

  gamma(2) += kappa*(modinv(2)-1);
  delta(4) += gamma(2) + modinv(2)*kappa;

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
