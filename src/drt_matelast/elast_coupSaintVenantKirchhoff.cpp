/*----------------------------------------------------------------------*/
/*!
\brief
This file contains the routines required for Saint-Venant-Kirchhoff material
The input line should read
  MAT 1 ELAST_CoupSVK YOUNG 1.044E7 NUE 0.3

\level 1

\maintainer Fabian Braeu
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
MAT::ELASTIC::PAR::CoupSVK::CoupSVK(Teuchos::RCP<MAT::PAR::Material> matdata) : Parameter(matdata)
{
  double c1 = matdata->GetDouble("YOUNG");
  double c2 = matdata->GetDouble("NUE");

  if (c2 <= 0.5 and c2 > -1.0)
  {
    lambda_ = (c2 == 0.5) ? 0.0 : c1 * c2 / ((1.0 + c2) * (1.0 - 2.0 * c2));
    mue_ = c1 / (2.0 * (1.0 + c2));  // shear modulus
  }
  else
    dserror("Poisson's ratio must be between -1.0 and 0.5!");
}

/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupSVK::CoupSVK(MAT::ELASTIC::PAR::CoupSVK* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupSVK::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  const double lambda = params_->lambda_;
  const double mue = params_->mue_;

  // strain energy: Psi = (1/4*mue+1/8*lambda)*I_1^2 - (0.75*lambda+0.5*mue)*I_1 - 0.5*mue*I_2 +
  // 9/8*lambda + 0.75*mue add to overall strain energy
  psi += (0.25 * mue + 0.125 * lambda) * prinv(0) * prinv(0) -
         (0.75 * lambda + 0.5 * mue) * prinv(0) - 0.5 * mue * prinv(1) + 1.125 * lambda +
         0.75 * mue;
}


/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupSVK::AddDerivativesPrincipal(LINALG::Matrix<3, 1>& dPI,
    LINALG::Matrix<6, 1>& ddPII, const LINALG::Matrix<3, 1>& prinv, const int eleGID)
{
  const double lambda = params_->lambda_;
  const double mue = params_->mue_;

  dPI(0) += (0.5 * mue + 0.25 * lambda) * prinv(0) - 0.75 * lambda - 0.5 * mue;
  dPI(1) -= 0.5 * mue;

  ddPII(0) += 0.5 * mue + 0.25 * lambda;


  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupSVK::AddThirdDerivativesPrincipalIso(
    LINALG::Matrix<10, 1>& dddPIII_iso, const LINALG::Matrix<3, 1>& prinv_iso, const int eleGID)
{
  // do nothing
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupSVK::AddCoupDerivVol(
    const double J, double* dPj1, double* dPj2, double* dPj3, double* dPj4)
{
  const double lambda = params_->lambda_;
  const double mu = params_->mue_;

  if (dPj1)
    *dPj1 += 12. * (mu / 4. + lambda / 8.) * pow(J, 1. / 3.) -
             2. * (3. / 4. * lambda + mu / 2.) * pow(J, -1. / 3.) - 2. * mu * pow(J, 1. / 3.);
  if (dPj2)
    *dPj2 += 4. * (mu / 4. + lambda / 8.) * pow(J, -2. / 3.) +
             2. / 3. * (3. / 4. * lambda + mu / 2.) * pow(J, -4. / 3.) -
             2. / 3. * mu * pow(J, -2. / 3.);
  if (dPj3)
    *dPj3 += -8. / 3. * (mu / 4. + lambda / 8.) * pow(J, -5. / 3.) -
             8. / 9. * (3. / 4. * lambda + mu / 2.) * pow(J, -7. / 3.) +
             4. / 9. * mu * pow(J, -5. / 3.);
  if (dPj4)
    *dPj4 += 40. / 9. * (mu / 4. + lambda / 8.) * pow(J, -8. / 3.) +
             56. / 27. * (3. / 4. * lambda + mu / 2.) * pow(J, -10. / 3.) -
             20. / 27. * mu * pow(J, -8. / 3.);
}

/*----------------------------------------------------------------------*/
