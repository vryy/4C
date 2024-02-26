/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a volumetric power law

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_volpow.hpp"

#include "baci_mat_par_material.hpp"

BACI_NAMESPACE_OPEN


MAT::ELASTIC::PAR::VolPow::VolPow(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), a_(*matdata->Get<double>("A")), expon_(*matdata->Get<double>("EXPON"))
{
}

MAT::ELASTIC::VolPow::VolPow(MAT::ELASTIC::PAR::VolPow* params) : params_(params) {}

void MAT::ELASTIC::VolPow::AddStrainEnergy(double& psi, const CORE::LINALG::Matrix<3, 1>& prinv,
    const CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  const double a = params_->a_;
  const double expon = params_->expon_;

  // strain energy: Psi = \frac{a}{expon-1} J^(-expon+1) + aJ
  // add to overall strain energy
  psi += a / (expon - 1.) * std::pow(modinv(2), -expon - 1.) + a * modinv(2);
}

void MAT::ELASTIC::VolPow::AddDerivativesModified(CORE::LINALG::Matrix<3, 1>& dPmodI,
    CORE::LINALG::Matrix<6, 1>& ddPmodII, const CORE::LINALG::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  const double a = params_->a_;
  const double expon = params_->expon_;

  dPmodI(2) += -a * (std::pow(modinv(2), -expon) - 1.);

  ddPmodII(2) += expon * a * std::pow(modinv(2), -(expon + 1.));
}

void MAT::ELASTIC::VolPow::Add3rdVolDeriv(
    const CORE::LINALG::Matrix<3, 1>& modinv, double& d3PsiVolDJ3)
{
  const double a = params_->a_;
  const double expon = params_->expon_;

  d3PsiVolDJ3 += -expon * a * (expon + 1.) * std::pow(modinv(2), -(expon + 2.));
}
BACI_NAMESPACE_CLOSE
