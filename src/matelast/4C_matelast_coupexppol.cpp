/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an isotropic exponential coupled material according to Weickenmeier_2014

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupexppol.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_material.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::PAR::CoupExpPol::CoupExpPol(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      a_(matdata->Get<double>("A")),
      b_(matdata->Get<double>("B")),
      c_(matdata->Get<double>("C"))
{
}

MAT::ELASTIC::CoupExpPol::CoupExpPol(MAT::ELASTIC::PAR::CoupExpPol* params) : params_(params) {}

void MAT::ELASTIC::CoupExpPol::AddStrainEnergy(double& psi, const CORE::LINALG::Matrix<3, 1>& prinv,
    const CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  const double a = params_->a_;
  const double b = params_->b_;
  const double c = params_->c_;

  // strain energy: Psi = a \exp[ b(I_1 - 3) - (2b + c)ln{J} + c(J-1) ] - a
  // add to overall strain energy
  psi += a * exp(b * (prinv(0) - 3.0) - (2.0 * b + c) * log(sqrt(prinv(2))) +
                 c * (sqrt(prinv(2)) - 1.0)) -
         a;
}

void MAT::ELASTIC::CoupExpPol::AddDerivativesPrincipal(CORE::LINALG::Matrix<3, 1>& dPI,
    CORE::LINALG::Matrix<6, 1>& ddPII, const CORE::LINALG::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  const double a = params_->a_;
  const double b = params_->b_;
  const double c = params_->c_;

  // ln of determinant of deformation gradient
  const double logdetf = std::log(std::sqrt(prinv(2)));

  // exponential function
  const double expfunc =
      std::exp(b * (prinv(0) - 3.0) - (2.0 * b + c) * logdetf + c * (std::sqrt(prinv(2)) - 1.0));

  dPI(0) += a * b * expfunc;
  dPI(2) += a * expfunc * (c / (2. * std::sqrt(prinv(2))) - (2. * b + c) / (2. * prinv(2)));

  ddPII(0) += a * b * b * expfunc;
  ddPII(2) += a * expfunc *
              ((0.5 * (2. * b + c) / (prinv(2) * prinv(2))) - (0.25 * c / std::pow(prinv(2), 1.5)) +
                  (-0.5 * (2. * b + c) / prinv(2) + 0.5 * c / std::sqrt(prinv(2))) *
                      (-0.5 * (2. * b + c) / prinv(2) + 0.5 * c / std::sqrt(prinv(2))));
  ddPII(4) += a * b * expfunc * (c / (2. * std::sqrt(prinv(2))) - (2. * b + c) / (2. * prinv(2)));
}
FOUR_C_NAMESPACE_CLOSE
