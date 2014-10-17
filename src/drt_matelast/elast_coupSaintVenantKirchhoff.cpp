/*----------------------------------------------------------------------*/
/*!
\file elast_coupSaintVenantKirchhoff.cpp
\brief


the input line should read either
  MAT 1 ELAST_CoupSVK YOUNG 1.044E7 NUE 0.3

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupSaintVenantKirchhoff.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupSVK::CoupSVK(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata)
{
  double c1 = matdata->GetDouble("YOUNG");
  double c2 = matdata->GetDouble("NUE");

  if (c2 <= 0.5 and c2 > -1.0)
  {
    lambda_ = (c2 == 0.5) ? 0.0 : c1*c2/((1.0+c2)*(1.0-2.0*c2));
    mue_ = c1/(2.0*(1.0+c2));  // shear modulus
  }
  else
    dserror("Poisson's ratio must be between -1.0 and 0.5!");
}

/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupSVK::CoupSVK(MAT::ELASTIC::PAR::CoupSVK* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupSVK::AddCoefficientsPrincipal(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{
  // gammas
  gamma(0) += .5*params_->lambda_*(prinv(0)-3.)-1.*params_->mue_;
  gamma(1) += params_->mue_;

  // deltas
  delta(0) += params_->lambda_;
  delta(7) += 2.*params_->mue_;

  return;
}


/*----------------------------------------------------------------------*/
