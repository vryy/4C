/*----------------------------------------------------------------------*/
/*!
\file elast_isoquad.cpp
\brief


the input line should read
  MAT 1 ELAST_IsoQuad C 1

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
#include "elast_isoquad.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoQuad::IsoQuad(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::IsoQuad::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::IsoQuad( this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoQuad::IsoQuad()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoQuad::IsoQuad(MAT::ELASTIC::PAR::IsoQuad* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoQuad::AddCoefficientsModified(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  const double c = params_ -> c_;

  //
  gamma(0) += 4*c*(modinv(0)-3);

  delta(0) += 8*c;

  return;
}


/*----------------------------------------------------------------------*/
