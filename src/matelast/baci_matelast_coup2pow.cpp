/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an isotropic general power-type material in terms of the second
Cauchy-Green invariant

\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_coup2pow.hpp"

#include "baci_mat_par_material.hpp"

BACI_NAMESPACE_OPEN

MAT::ELASTIC::PAR::Coup2Pow::Coup2Pow(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), c_(matdata->GetDouble("C")), d_(matdata->GetInt("D"))
{
}

MAT::ELASTIC::Coup2Pow::Coup2Pow(MAT::ELASTIC::PAR::Coup2Pow* params) : params_(params) {}

void MAT::ELASTIC::Coup2Pow::AddStrainEnergy(double& psi, const CORE::LINALG::Matrix<3, 1>& prinv,
    const CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const int d = params_->d_;

  // strain energy: Psi = C (II_{\boldsymbol{C}}-3)^D
  // add to overall strain energy
  psi += c * pow((prinv(1) - 3.), d);
}

void MAT::ELASTIC::Coup2Pow::AddDerivativesPrincipal(CORE::LINALG::Matrix<3, 1>& dPI,
    CORE::LINALG::Matrix<6, 1>& ddPII, const CORE::LINALG::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  const double c = params_->c_;
  const int d = params_->d_;

  // If d<2 the material model is not stress free in the reference configuration
  if (d < 2)
    dserror(
        "The Elast_Coup2Pow - material only works for positive integer exponents, which are larger "
        "than two.");

  dPI(1) += c * d * pow((prinv(1) - 3.), d - 1.);

  if (d == 2)
    ddPII(1) += (c * d * d - c * d);
  else
    ddPII(1) += (c * d * d - c * d) * pow((prinv(1) - 3.), d - 2.);
}
BACI_NAMESPACE_CLOSE
