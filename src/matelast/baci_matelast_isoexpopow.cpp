/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the isochoric contribution of an isotropic exponential material

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_isoexpopow.hpp"

#include "baci_mat_par_material.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::PAR::IsoExpoPow::IsoExpoPow(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      k1_(*matdata->Get<double>("K1")),
      k2_(*matdata->Get<double>("K2")),
      d_(*matdata->Get<int>("C"))
{
}

MAT::ELASTIC::IsoExpoPow::IsoExpoPow(MAT::ELASTIC::PAR::IsoExpoPow* params) : params_(params) {}

void MAT::ELASTIC::IsoExpoPow::AddStrainEnergy(double& psi, const CORE::LINALG::Matrix<3, 1>& prinv,
    const CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;
  const int d = params_->d_;

  // strain energy: Psi = \frac{k_1}{2k_2} (e^{k_2 (\overline{I}_{\boldsymbol{C}}-3)^d}-1).
  // add to overall strain energy
  if (k2 != 0) psi += k1 / (2. * k2) * (exp(k2 * pow(modinv(0) - 3., d)) - 1.);
}

void MAT::ELASTIC::IsoExpoPow::AddDerivativesModified(CORE::LINALG::Matrix<3, 1>& dPmodI,
    CORE::LINALG::Matrix<6, 1>& ddPmodII, const CORE::LINALG::Matrix<3, 1>& modinv, const int gp,
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
}
FOUR_C_NAMESPACE_CLOSE
