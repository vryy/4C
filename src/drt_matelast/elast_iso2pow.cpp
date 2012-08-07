/*----------------------------------------------------------------------*/
/*!
\file elast_iso2pow.cpp
\brief


the input line should read
  MAT 1 ELAST_Iso2Pow C 1 D 1

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_iso2pow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Iso2Pow::Iso2Pow(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  d_(matdata->GetInt("D"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::Iso2Pow::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::Iso2Pow( this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Iso2Pow::Iso2Pow()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Iso2Pow::Iso2Pow(MAT::ELASTIC::PAR::Iso2Pow* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Iso2Pow::AddCoefficientsModified(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{

  const double c = params_ -> c_;
  const    int d = params_ -> d_;

  if (d<=0)
    dserror("The Elast_Iso2Pow - material only works for positive integer exponents larger than one.");


  gamma(0) += 2.*modinv(0)*c*d*pow((modinv(1)-3),d-1);
  gamma(1) -= 2.*c*d*pow((modinv(1)-3),d-1);

  delta(0) += 4*c*d*pow((modinv(1)-3),d-1);
  delta(3) -= 4*c*d*pow((modinv(1)-3),d-1);

  if (d>=2)
  {
    delta(0) += 4.*modinv(0)*modinv(0)*c*d*(d-1)*pow((modinv(1)-3),d-2);
    delta(1) -= 4.*modinv(0)*c*d*(d-1)*pow((modinv(1)-3),d-2);
    delta(2) += 4*c*d*(d-1)*pow((modinv(1)-3),d-2);
  }

  else
    delta(0) += 0;


  return;
}


/*----------------------------------------------------------------------*/
