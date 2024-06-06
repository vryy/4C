/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the volumetic contribution suggested by Ogden, see "Doll, S. and
Schweizerhof, K. On the Development of Volumetric Strain Energy Functions, Journal of Applied
Mechanics, 2000"

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_vologden.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::VolOgden::VolOgden(const Teuchos::RCP<Core::Mat::PAR::Material>& matdata)
    : Parameter(matdata), kappa_(matdata->Get<double>("KAPPA")), beta_(matdata->Get<double>("BETA"))
{
}

Mat::Elastic::VolOgden::VolOgden(Mat::Elastic::PAR::VolOgden* params) : params_(params) {}

void Mat::Elastic::VolOgden::AddStrainEnergy(double& psi, const Core::LinAlg::Matrix<3, 1>& prinv,
    const Core::LinAlg::Matrix<3, 1>& modinv, const Core::LinAlg::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  const double kappa = params_->kappa_;
  const double beta = params_->beta_;

  // strain energy: Psi = \frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)
  // add to overall strain energy
  if (beta != 0)
    psi += kappa / (beta * beta) * (beta * log(modinv(2)) + pow(modinv(2), -beta) - 1.);
  else
    psi += kappa / 2. * pow(std::log(modinv(2)), 2.);
}

void Mat::Elastic::VolOgden::add_derivatives_modified(Core::LinAlg::Matrix<3, 1>& dPmodI,
    Core::LinAlg::Matrix<6, 1>& ddPmodII, const Core::LinAlg::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  const double kappa = params_->kappa_;
  const double beta = params_->beta_;

  if (beta != 0.)
  {
    dPmodI(2) += kappa / (beta * beta) * (beta / modinv(2) - beta * pow(modinv(2), -beta - 1.));
    ddPmodII(2) +=
        kappa / (beta * beta) *
        (-beta / (modinv(2) * modinv(2)) + beta * (beta + 1.) * pow(modinv(2), -beta - 2.));
  }
  else
  {
    dPmodI(2) += kappa * std::log(modinv(2)) / modinv(2);
    ddPmodII(2) += kappa * (1. - std::log(modinv(2))) / (modinv(2) * modinv(2));
  }
}

void Mat::Elastic::VolOgden::Add3rdVolDeriv(
    const Core::LinAlg::Matrix<3, 1>& modinv, double& d3PsiVolDJ3)
{
  const double kappa = params_->kappa_;
  const double beta = params_->beta_;
  const double J = modinv(2);
  if (beta != 0.)
    d3PsiVolDJ3 +=
        kappa / (beta * J * J * J) * (2. - pow(J, -beta) * (beta * beta + 3. * beta + 2.));
  else
    d3PsiVolDJ3 += kappa * (-3. + 2. * std::log(modinv(2))) / (modinv(2) * modinv(2) * modinv(2));
}
FOUR_C_NAMESPACE_CLOSE
