/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an isotropic general power-type material in terms of the first Cauchy-Green
invariant

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_coup1pow.hpp"

#include "baci_mat_par_material.hpp"

FOUR_C_NAMESPACE_OPEN

MAT::ELASTIC::PAR::Coup1Pow::Coup1Pow(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), c_(*matdata->Get<double>("C")), d_(*matdata->Get<int>("D"))
{
}

MAT::ELASTIC::Coup1Pow::Coup1Pow(MAT::ELASTIC::PAR::Coup1Pow* params) : params_(params) {}

void MAT::ELASTIC::Coup1Pow::AddStrainEnergy(double& psi, const CORE::LINALG::Matrix<3, 1>& prinv,
    const CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const int d = params_->d_;

  // strain energy: Psi = C (I_{\boldsymbol{C}}-3)^D
  // add to overall strain energy
  psi += c * pow((prinv(0) - 3.), d);
}

void MAT::ELASTIC::Coup1Pow::AddDerivativesPrincipal(CORE::LINALG::Matrix<3, 1>& dPI,
    CORE::LINALG::Matrix<6, 1>& ddPII, const CORE::LINALG::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  const double c = params_->c_;
  const double d = params_->d_;

  // If d<2 the material model is not stress free in the reference configuration
  if (d < 2)
    dserror(
        "The Elast_Coup1Pow - material only works for positive integer exponents, which are larger "
        "than two.");

  dPI(0) += c * d * pow((prinv(0) - 3.), d - 1.);

  if (d == 2)
    ddPII(0) += (c * d * d - c * d);
  else
    ddPII(0) += (c * d * d - c * d) * pow((prinv(0) - 3.), d - 2.);
}
FOUR_C_NAMESPACE_CLOSE
