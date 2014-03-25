/*----------------------------------------------------------------------*/
/*!
\file elast_coupblatzko.cpp
\brief


the input line should read
  MAT 1 ELAST_CoupBlatzKo MUE 1.044E7 NUE 0.3 F 0.5

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupblatzko.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupBlatzKo::CoupBlatzKo(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  mue_(matdata->GetDouble("MUE")),
  nue_(matdata->GetDouble("NUE")),
  f_(matdata->GetDouble("F"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupBlatzKo::CoupBlatzKo(MAT::ELASTIC::PAR::CoupBlatzKo* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupBlatzKo::AddCoefficientsPrincipal(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{

  // material parameters for isochoric part
  const double mue  = params_->mue_;  // Shear modulus
  const double nue  = params_->nue_;  // Poisson's ratio
  const double f    = params_->f_;    // interpolation parameter

  // introducing xeta for consitence with Holzapfel and simplification
  const double beta = nue/(1. - 2.*nue);
  const double xeta = 2.*mue*(1.-f)/prinv(2);

  // gammas
  gamma(0) += mue*f+xeta*prinv(0)/2.;
  gamma(1) += -xeta/2.;
  gamma(2) += -mue*f*pow(prinv(2), -beta) - xeta*(prinv(1)-pow(prinv(2), beta+1))/2.;


  // deltas
  delta(0) += xeta;
  delta(2) += -xeta*prinv(0);
  delta(4) += xeta;
  delta(5) += 2*mue*f*beta*pow(prinv(2), -beta) + xeta*(prinv(1) + beta*pow(prinv(2), (beta+1)));
  delta(6) += 2*mue*f*pow(prinv(2), -beta) + xeta*(prinv(1) - pow(prinv(2), (beta+1)));
  delta(7) += -xeta;

  return;
}


/*----------------------------------------------------------------------*/
