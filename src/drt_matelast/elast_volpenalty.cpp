/*----------------------------------------------------------------------*/
/*!
\file elast_volpenalty.cpp
\brief


the input line should read
  MAT 1 ELAST_VolPenalty EPSILON 1. GAMMA 1.

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
#include "elast_volpenalty.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::VolPenalty::VolPenalty(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  eps_(matdata->GetDouble("EPSILON")),
  gam_(matdata->GetDouble("GAMMA"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VolPenalty::VolPenalty(MAT::ELASTIC::PAR::VolPenalty* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolPenalty::AddCoefficientsModified(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  const double eps = params_ -> eps_;
  const double gam = params_ -> gam_;

  gamma(2) += eps*gam * (pow(modinv(2), gam)-pow(modinv(2), - gam)) / modinv(2);  ;
  delta(4) +=  eps*gam*gam * (pow(modinv(2), gam)+pow(modinv(2), - gam)) / modinv(2);  ;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolPenalty::Add3rdVolDeriv(const LINALG::Matrix<3,1>& modinv, double& d3PsiVolDJ3)
{
  const double eps = params_ -> eps_;
  const double gam = params_ -> gam_;
  const double J=modinv(2);
  d3PsiVolDJ3 += eps*(   gam *( gam-1.)*( gam-2.)*pow(J, gam-3.) +
                       (-gam)*(-gam-1.)*(-gam-2.)*pow(J,-gam-3.)
                     );
  return;
}


/*----------------------------------------------------------------------*/
