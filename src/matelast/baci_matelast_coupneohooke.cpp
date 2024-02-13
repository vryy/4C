/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a coupled Neo Hookean material


\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_coupneohooke.hpp"

#include "baci_mat_par_material.hpp"

#include <limits>

BACI_NAMESPACE_OPEN


MAT::ELASTIC::PAR::CoupNeoHooke::CoupNeoHooke(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), youngs_(matdata->GetDouble("YOUNG")), nue_(matdata->GetDouble("NUE"))
{
  // Material Constants c and beta
  c_ = youngs_ / (4.0 * (1.0 + nue_));
  beta_ = nue_ / (1.0 - 2.0 * nue_);
}

MAT::ELASTIC::CoupNeoHooke::CoupNeoHooke(MAT::ELASTIC::PAR::CoupNeoHooke* params) : params_(params)
{
}

void MAT::ELASTIC::CoupNeoHooke::AddShearMod(
    bool& haveshearmod,  ///< non-zero shear modulus was added
    double& shearmod     ///< variable to add upon
) const
{
  haveshearmod = true;

  shearmod += 2 * params_->c_;
}

void MAT::ELASTIC::CoupNeoHooke::AddStrainEnergy(double& psi,
    const CORE::LINALG::Matrix<3, 1>& prinv, const CORE::LINALG::Matrix<3, 1>& modinv,
    const CORE::LINALG::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const double beta = params_->beta_;

  // strain energy: psi = c / beta * (I3^{-beta} - 1) + c * (I1 - 3)
  double psiadd = c * (prinv(0) - 3.);
  if (beta != 0)  // take care of possible division by zero in case or Poisson's ratio nu = 0.0
    psiadd += (c / beta) * (std::exp(std::log(prinv(2)) * (-beta)) - 1.);
  else
    psiadd -= c * std::log(prinv(2));

  // add to overall strain energy
  psi += psiadd;
}

void MAT::ELASTIC::CoupNeoHooke::AddDerivativesPrincipal(CORE::LINALG::Matrix<3, 1>& dPI,
    CORE::LINALG::Matrix<6, 1>& ddPII, const CORE::LINALG::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  const double beta = params_->beta_;
  const double c = params_->c_;

  dPI(0) += c;
  // computing exp(log(a)*b) is faster than pow(a,b)
  if (prinv(2) > 0)
  {
    const double prinv2_to_beta_m1 = std::exp(std::log(prinv(2)) * (-beta - 1.));
    dPI(2) -= c * prinv2_to_beta_m1;
    ddPII(2) += c * (beta + 1.) * prinv2_to_beta_m1 / prinv(2);
  }
  else
    dPI(2) = ddPII(2) = std::numeric_limits<double>::quiet_NaN();
}

void MAT::ELASTIC::CoupNeoHooke::AddThirdDerivativesPrincipalIso(
    CORE::LINALG::Matrix<10, 1>& dddPIII_iso, const CORE::LINALG::Matrix<3, 1>& prinv_iso,
    const int gp, const int eleGID)
{
  const double beta = params_->beta_;
  const double c = params_->c_;

  dddPIII_iso(2) -= c * (beta + 1.0) * (beta + 2.0) * std::pow(prinv_iso(2), -beta - 3.0);
}

void MAT::ELASTIC::CoupNeoHooke::AddCoupDerivVol(
    const double J, double* dPj1, double* dPj2, double* dPj3, double* dPj4)
{
  const double beta = params_->beta_;
  const double c = params_->c_;

  if (J < 0.) dserror("negative jacobian determinant");

  if (dPj1) *dPj1 += 2. * c * pow(J, -1. / 3.) - 2. * c * pow(J * J, -beta) / J;
  if (dPj2)
    *dPj2 += -2. / 3. * c * pow(J, -4. / 3.) + 4. * c * pow(J * J, -beta) * beta * pow(J, -2.) +
             2. * c * pow(J * J, -beta) * pow(J, -2.);
  if (dPj3)
    *dPj3 += 0.8e1 / 0.9e1 * c * pow(J, -0.7e1 / 3.) -
             0.8e1 * c * pow(J * J, -beta) * beta * beta * pow(J, -3.) -
             0.12e2 * c * pow(J * J, -beta) * beta * pow(J, -3.) -
             4. * c * pow(J * J, -beta) * pow(J, -3.);
  if (dPj4)
    *dPj4 += -56. / 27. * c * pow(J, -10. / 3.) +
             16. * c * pow(J * J, -beta) * pow(beta, 3.) * pow(J, -4.) +
             48. * c * pow(J * J, -beta) * beta * beta * pow(J, -4.) +
             44. * c * pow(J * J, -beta) * beta * pow(J, -4.) +
             12. * c * pow(J * J, -beta) * pow(J, -4.);
}
BACI_NAMESPACE_CLOSE
