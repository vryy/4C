/*----------------------------------------------------------------------*/
/*!
\file elast_coup2pow.cpp
\brief


the input line should read
  MAT 1 ELAST_Coup2Pow C 1 D 1

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
#include "elast_coup2pow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Coup2Pow::Coup2Pow(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  d_(matdata->GetInt("D"))
{
}


/*----------------------------------------------------------------------*
 *            Constructor Material Class                               *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Coup2Pow::Coup2Pow(MAT::ELASTIC::PAR::Coup2Pow* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup2Pow::AddCoefficientsPrincipal(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{

  const double c = params_ -> c_;
  const    int d = params_ -> d_;


  // If d<2 the material model is not stress free in the reference configuration
  if (d<2)
    dserror("The Elast_Coup2Pow - material only works for positive integer exponents, which are larger than two.");


  gamma(0) += 2.*prinv(0)*c*d*pow((prinv(1)-3),d-1);
  gamma(1) -= 2.*c*d*pow((prinv(1)-3),d-1);

  delta(0) += 4*c*d*pow((prinv(1)-3),d-1)+4.*prinv(0)*prinv(0)*c*d*(d-1)*pow((prinv(1)-3),d-2);
  delta(1) -= 4.*prinv(0)*c*d*(d-1)*pow((prinv(1)-3),d-2);
  delta(3) += 4*c*d*(d-1)*pow((prinv(1)-3),d-2);
  delta(7) -= 4*c*d*pow((prinv(1)-3),d-1);


  return;
}


/*----------------------------------------------------------------------*/
