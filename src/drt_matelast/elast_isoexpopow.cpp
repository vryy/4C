/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required to calculate the isochoric contribution
of an exponential material.
The input line should read
MAT 1 ELAST_IsoExpoPow K1 5000. K2 5. C 1.
C = Exponent D

\level 1

*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isoexpopow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoExpoPow::IsoExpoPow(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      k1_(matdata->GetDouble("K1")),
      k2_(matdata->GetDouble("K2")),
      d_(matdata->GetInt("C"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoExpoPow::IsoExpoPow(MAT::ELASTIC::PAR::IsoExpoPow* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoExpoPow::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int gp,
    const int eleGID)
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;
  const int d = params_->d_;

  // strain energy: Psi = \frac{k_1}{2k_2} (e^{k_2 (\overline{I}_{\boldsymbol{C}}-3)^d}-1).
  // add to overall strain energy
  if (k2 != 0) psi += k1 / (2. * k2) * (exp(k2 * pow(modinv(0) - 3., d)) - 1.);
}


/*----------------------------------------------------------------------
 *                                                      birzle 11/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoExpoPow::AddDerivativesModified(LINALG::Matrix<3, 1>& dPmodI,
    LINALG::Matrix<6, 1>& ddPmodII, const LINALG::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;
  const int d = params_->d_;
  const double expf = std::exp(k2 * std::pow(modinv(0) - 3., (double)d));

  if (d < 1)
    dserror(
        "The Elast_IsoExpoPow - material only works for positive integer exponents (C) larger than "
        "one.");

  if (d == 1)
    dPmodI(0) += k1 * d * 0.5 * expf;
  else
    dPmodI(0) += k1 * d * 0.5 * expf * pow(modinv(0) - 3., d - 1.);

  if (d == 1)
    ddPmodII(0) += k1 * d * d * k2 * 0.5 * expf;
  else if (d == 2)
    ddPmodII(0) += k1 * d * d * k2 * 0.5 * expf * pow((modinv(0) - 3.), 2. * (d - 1.)) +
                   k1 * d * 0.5 * expf * (d - 1.);
  else
    ddPmodII(0) += k1 * d * d * k2 * 0.5 * expf * pow((modinv(0) - 3.), 2. * (d - 1.)) +
                   k1 * d * 0.5 * expf * (d - 1.) * pow((modinv(0) - 3.), (d - 2.));

  return;
}
