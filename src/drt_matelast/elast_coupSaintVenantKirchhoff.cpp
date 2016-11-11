/*----------------------------------------------------------------------*/
/*!
\file elast_coupSaintVenantKirchhoff.cpp
\brief the input line should read either
  MAT 1 ELAST_CoupSVK YOUNG 1.044E7 NUE 0.3

\level 1

<pre>
\maintainer Alexander Seitz
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
void MAT::ELASTIC::CoupSVK::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1> glstrain,
    const int eleGID)
{
  const double lambda  = params_->lambda_;
  const double mue     = params_->mue_;

  // strain energy: Psi = (1/4*mue+1/8*lambda)*I_1^2 - (0.75*lambda+0.5*mue)*I_1 - 0.5*mue*I_2 + 9/8*lambda + 0.75*mue
  // add to overall strain energy
  psi += (0.25*mue+0.125*lambda)*prinv(0)*prinv(0) - (0.75*lambda+0.5*mue)*prinv(0) - 0.5*mue*prinv(1) + 1.125*lambda + 0.75*mue;
}


/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupSVK::AddDerivativesPrincipal(
    LINALG::Matrix<3,1>& dPI,
    LINALG::Matrix<6,1>& ddPII,
    const LINALG::Matrix<3,1>& prinv,
    const int eleGID
)
{
  const double lambda  = params_->lambda_;
  const double mue     = params_->mue_;

  dPI(0) += (0.5*mue+0.25*lambda)*prinv(0)-0.75*lambda-0.5*mue;
  dPI(1) -= 0.5*mue;

  ddPII(0) += 0.5*mue+0.25*lambda;


  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupSVK::AddThirdDerivativesPrincipalIso(LINALG::Matrix<10,1>& dddPIII_iso,
                                                                 const LINALG::Matrix<3,1>& prinv_iso,
                                                                 const int eleGID)
{
// do nothing
}


/*----------------------------------------------------------------------*/
