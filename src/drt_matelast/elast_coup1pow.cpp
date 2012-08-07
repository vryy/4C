/*----------------------------------------------------------------------*/
/*!
\file elast_coup1pow.cpp
\brief


the input line should read
  MAT 1 ELAST_Coup1Pow C 1 D 1

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
#include "elast_coup1pow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Coup1Pow::Coup1Pow(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  d_(matdata->GetInt("D"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::Coup1Pow::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::Coup1Pow( this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Coup1Pow::Coup1Pow()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Coup1Pow::Coup1Pow(MAT::ELASTIC::PAR::Coup1Pow* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup1Pow::AddCoefficientsPrincipal(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{

  const double c = params_ -> c_;
  const    int d = params_ -> d_;


  // If d<2 the material model is not stress free in the reference configuration
  if (d<2)
    dserror("The Elast_Coup1Pow - material only works for positive integer exponents, which are larger than two.");


  gamma(0) += 2.*c*d*pow((prinv(0)-3),d-1);

  delta(0) += 4.*c*d*(d-1)*pow((prinv(0)-3),d-2);


  return;
}


/*----------------------------------------------------------------------*/
