/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the volumetric SussmanBathe material according to "Doll, S. and
Schweizerhof, K. On the Development of Volumetric Strain Energy Functions Journal of Applied
Mechanics, 2000"


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_volsussmanbathe.hpp"

#include "4C_material_input_base.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::PAR::VolSussmanBathe::VolSussmanBathe(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : Parameter(matdata), kappa_(matdata->Get<double>("KAPPA"))
{
}

MAT::ELASTIC::VolSussmanBathe::VolSussmanBathe(MAT::ELASTIC::PAR::VolSussmanBathe* params)
    : params_(params)
{
}

void MAT::ELASTIC::VolSussmanBathe::AddStrainEnergy(double& psi,
    const CORE::LINALG::Matrix<3, 1>& prinv, const CORE::LINALG::Matrix<3, 1>& modinv,
    const CORE::LINALG::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  const double kappa = params_->kappa_;

  // strain energy: Psi = \frac \kappa 2 (J-1)^2
  // add to overall strain energy
  psi += kappa * 0.5 * (modinv(2) - 1.) * (modinv(2) - 1.);
}

void MAT::ELASTIC::VolSussmanBathe::AddDerivativesModified(CORE::LINALG::Matrix<3, 1>& dPmodI,
    CORE::LINALG::Matrix<6, 1>& ddPmodII, const CORE::LINALG::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  const double kappa = params_->kappa_;

  dPmodI(2) += kappa * (modinv(2) - 1.);

  ddPmodII(2) += kappa;
}

void MAT::ELASTIC::VolSussmanBathe::Add3rdVolDeriv(
    const CORE::LINALG::Matrix<3, 1>& modinv, double& d3PsiVolDJ3)
{
  d3PsiVolDJ3 += 0.;
}
FOUR_C_NAMESPACE_CLOSE
