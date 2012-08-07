/*----------------------------------------------------------------------*/
/*!
\file elast_iso1pow.cpp
\brief


the input line should read
  MAT 1 ELAST_Iso1Pow C 1 D 1

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
#include "elast_iso1pow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Iso1Pow::Iso1Pow(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  d_(matdata->GetInt("D"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::Iso1Pow::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::Iso1Pow( this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Iso1Pow::Iso1Pow()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Iso1Pow::Iso1Pow(MAT::ELASTIC::PAR::Iso1Pow* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Iso1Pow::AddCoefficientsModified(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{

  const double c = params_ -> c_;
  const    int d = params_ -> d_;

  if (d<=0)
    dserror("The Elast_Iso1Pow - material only works for positive integer exponents larger than one.");


  gamma(0) += 2.*c*d*pow((modinv(0)-3),d-1);

  if (d>=2)
    delta(0) += 4.*c*d*(d-1)*pow((modinv(0)-3),d-2);
  else
    delta(0) += 0;


  return;
}


/*----------------------------------------------------------------------*/
