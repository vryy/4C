/*----------------------------------------------------------------------*/
/*!
\file elast_isoyeoh.cpp

\brief
This file contains the routines required to calculate the isochoric contribution
of a Yeoh-type material.
The input line should read
  MAT 1 ELAST_IsoYeoh C1 100 C2 0 C3 200

\level 1

<pre>
\maintainer Fabian Braeu
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isoyeoh.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoYeoh::IsoYeoh(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      c1_(matdata->GetDouble("C1")),
      c2_(matdata->GetDouble("C2")),
      c3_(matdata->GetDouble("C3"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoYeoh::IsoYeoh(MAT::ELASTIC::PAR::IsoYeoh* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoYeoh::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;
  const double c3 = params_->c3_;

  // strain energy: Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + C2
  // (\overline{I}_{\boldsymbol{C}}-3)^2 + C3 (\overline{I}_{\boldsymbol{C}}-3)^3. add to overall
  // strain energy
  psi += c1 * (modinv(0) - 3.) + c2 * (modinv(0) - 3.) * (modinv(0) - 3.) +
         c3 * (modinv(0) - 3.) * (modinv(0) - 3.) * (modinv(0) - 3.);
}


/*----------------------------------------------------------------------
 *                                                      birzle 12/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoYeoh::AddDerivativesModified(LINALG::Matrix<3, 1>& dPmodI,
    LINALG::Matrix<6, 1>& ddPmodII, const LINALG::Matrix<3, 1>& modinv, const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;
  const double c3 = params_->c3_;

  dPmodI(0) += c1 + 2. * c2 * (modinv(0) - 3.) + 3. * c3 * (modinv(0) - 3.) * (modinv(0) - 3.);
  ddPmodII(0) += 2. * c2 + 6. * c3 * (modinv(0) - 3.);

  return;
}


/*----------------------------------------------------------------------*/
