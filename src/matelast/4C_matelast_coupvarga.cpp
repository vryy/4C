/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the isotropic Varga material
\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupvarga.hpp"

#include "4C_mat_par_material.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::PAR::CoupVarga::CoupVarga(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), mue_(*matdata->Get<double>("MUE")), beta_(*matdata->Get<double>("BETA"))
{
}

MAT::ELASTIC::CoupVarga::CoupVarga(MAT::ELASTIC::PAR::CoupVarga* params) : params_(params) {}

void MAT::ELASTIC::CoupVarga::AddShearMod(bool& haveshearmod,  ///< non-zero shear modulus was added
    double& shearmod                                           ///< variable to add upon
) const
{
  // indeed, a shear modulus is provided
  haveshearmod = haveshearmod or true;

  // material parameters for isochoric part
  shearmod += params_->mue_;
}

void MAT::ELASTIC::CoupVarga::AddCoefficientsStretchesPrincipal(
    CORE::LINALG::Matrix<3, 1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
    CORE::LINALG::Matrix<6, 1>&
        delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23, delta_31]
    const CORE::LINALG::Matrix<3, 1>&
        prstr  ///< principal stretches, [lambda_1, lambda_2, lambda_3]
)
{
  // parameters
  const double alpha = 2.0 * params_->mue_ - params_->beta_;
  const double beta = params_->beta_;

  // first derivatives
  // \frac{\partial Psi}{\partial \lambda_1}
  gamma(0) += alpha - beta / (prstr(0) * prstr(0));
  // \frac{\partial Psi}{\partial \lambda_2}
  gamma(1) += alpha - beta / (prstr(1) * prstr(1));
  // \frac{\partial Psi}{\partial \lambda_3}
  gamma(2) += alpha - beta / (prstr(2) * prstr(2));

  // second derivatives
  // \frac{\partial^2 Psi}{\partial\lambda_1^2}
  delta(0) += 2.0 * beta / (prstr(0) * prstr(0) * prstr(0));
  // \frac{\partial^2 Psi}{\partial\lambda_2^2}
  delta(1) += 2.0 * beta / (prstr(1) * prstr(1) * prstr(1));
  // \frac{\partial^2 Psi}{\partial\lambda_3^2}
  delta(2) += 2.0 * beta / (prstr(2) * prstr(2) * prstr(2));
  // \frac{\partial^2 Psi}{\partial\lambda_1 \partial\lambda_2}
  delta(3) += 0.0;
  // \frac{\partial^2 Psi}{\partial\lambda_2 \partial\lambda_3}
  delta(4) += 0.0;
  // \frac{\partial^2 Psi}{\partial\lambda_3 \partial\lambda_1}
  delta(5) += 0.0;
}
FOUR_C_NAMESPACE_CLOSE
