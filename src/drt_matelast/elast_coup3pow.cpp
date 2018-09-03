/*----------------------------------------------------------------------*/
/*!
\file elast_coup3pow.cpp
\brief the input line should read
  MAT 1 ELAST_Coup3Pow C 1 D 1

\level 2

<pre>
\maintainer Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coup3pow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Coup3Pow::Coup3Pow(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), c_(matdata->GetDouble("C")), d_(matdata->GetInt("D"))
{
}

/*----------------------------------------------------------------------*
 *            Constructor Material Class                                *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Coup3Pow::Coup3Pow(MAT::ELASTIC::PAR::Coup3Pow* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup3Pow::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1> glstrain, const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const int d = params_->d_;

  // add to overall strain energy
  psi += c * pow((pow(prinv(2), 1. / 3.) - 1.), d);
}


/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup3Pow::AddDerivativesPrincipal(LINALG::Matrix<3, 1>& dPI,
    LINALG::Matrix<6, 1>& ddPII, const LINALG::Matrix<3, 1>& prinv, const int eleGID)
{
  const double c = params_->c_;
  const int d = params_->d_;

  // If d<2
  if (d < 2)
    dserror(
        "The Elast_Coup3Pow - material only works for positive integer exponents, which are larger "
        "than two.");

  dPI(2) += 1. / 3. * c * d * pow(prinv(2), -2. / 3.) * pow((pow(prinv(2), 1. / 3.) - 1.), d - 1.);

  if (d == 2)
    ddPII(2) +=
        -2. / 9. * c * d * pow(prinv(2), -5. / 3.) * pow((pow(prinv(2), 1. / 3.) - 1.), d - 1.) +
        1. / 9. * c * d * (d - 1.) * pow(prinv(2), -4. / 3.);
  else
    ddPII(2) +=
        -2. / 9. * c * d * pow(prinv(2), -5. / 3.) * pow((pow(prinv(2), 1. / 3.) - 1.), d - 1.) +
        1. / 9. * c * d * (d - 1.) * pow(prinv(2), -4. / 3.) *
            pow((pow(prinv(2), 1. / 3.) - 1.), d - 2.);

  return;
}

/*----------------------------------------------------------------------*/
