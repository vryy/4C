/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the volumetric SussmanBathe material according to "Doll, S. and
Schweizerhof, K. On the Development of Volumetric Strain Energy Functions Journal of Applied
Mechanics, 2000"


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_volsussmanbathe.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::VolSussmanBathe::VolSussmanBathe(
    const Teuchos::RCP<Core::Mat::PAR::Material>& matdata)
    : Parameter(matdata), kappa_(matdata->Get<double>("KAPPA"))
{
}

Mat::Elastic::VolSussmanBathe::VolSussmanBathe(Mat::Elastic::PAR::VolSussmanBathe* params)
    : params_(params)
{
}

void Mat::Elastic::VolSussmanBathe::AddStrainEnergy(double& psi,
    const Core::LinAlg::Matrix<3, 1>& prinv, const Core::LinAlg::Matrix<3, 1>& modinv,
    const Core::LinAlg::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  const double kappa = params_->kappa_;

  // strain energy: Psi = \frac \kappa 2 (J-1)^2
  // add to overall strain energy
  psi += kappa * 0.5 * (modinv(2) - 1.) * (modinv(2) - 1.);
}

void Mat::Elastic::VolSussmanBathe::add_derivatives_modified(Core::LinAlg::Matrix<3, 1>& dPmodI,
    Core::LinAlg::Matrix<6, 1>& ddPmodII, const Core::LinAlg::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  const double kappa = params_->kappa_;

  dPmodI(2) += kappa * (modinv(2) - 1.);

  ddPmodII(2) += kappa;
}

void Mat::Elastic::VolSussmanBathe::Add3rdVolDeriv(
    const Core::LinAlg::Matrix<3, 1>& modinv, double& d3PsiVolDJ3)
{
  d3PsiVolDJ3 += 0.;
}
FOUR_C_NAMESPACE_CLOSE
