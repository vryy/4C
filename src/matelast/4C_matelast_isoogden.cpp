/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the isochoric one-term Ogden material.

\level 3
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_isoogden.hpp"

#include "4C_material_input_base.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::PAR::IsoOgden::IsoOgden(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata), mue_(matdata->Get<double>("MUE")), alpha_(matdata->Get<double>("ALPHA"))
{
}

MAT::ELASTIC::IsoOgden::IsoOgden(MAT::ELASTIC::PAR::IsoOgden* params) : params_(params) {}

void MAT::ELASTIC::IsoOgden::AddCoefficientsStretchesModified(CORE::LINALG::Matrix<3, 1>& modgamma,
    CORE::LINALG::Matrix<6, 1>& moddelta, const CORE::LINALG::Matrix<3, 1>& modstr)
{
  // get parameters
  const double& mue = params_->mue_;
  const double& alpha = params_->alpha_;

  // first derivatives \frac{\partial Psi}{\partial \bar{\lambda}_\i}
  modgamma(0) += 2 * mue / alpha * std::pow(modstr(0), alpha - 1);  // ,0
  modgamma(1) += 2 * mue / alpha * std::pow(modstr(1), alpha - 1);  // ,1
  modgamma(2) += 2 * mue / alpha * std::pow(modstr(2), alpha - 1);  // ,2

  // second derivatives \frac{\partial^2 Psi}{\partial \bar{\lambda}_\i \partial \bar {\lambda}_\j}
  moddelta(0) += 2 * mue / alpha * (alpha - 1) * std::pow(modstr(0), alpha - 2);  // ,00
  moddelta(1) += 2 * mue / alpha * (alpha - 1) * std::pow(modstr(1), alpha - 2);  // ,11
  moddelta(2) += 2 * mue / alpha * (alpha - 1) * std::pow(modstr(2), alpha - 2);  // ,22
  moddelta(3) += 0;                                                               // ,01
  moddelta(4) += 0;                                                               // ,12
  moddelta(5) += 0;                                                               // ,20
}

FOUR_C_NAMESPACE_CLOSE
