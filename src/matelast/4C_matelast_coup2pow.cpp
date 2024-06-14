/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an isotropic general power-type material in terms of the second
Cauchy-Green invariant

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coup2pow.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::Elastic::PAR::Coup2Pow::Coup2Pow(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      c_(matdata.parameters.Get<double>("C")),
      d_(matdata.parameters.Get<int>("D"))
{
}

Mat::Elastic::Coup2Pow::Coup2Pow(Mat::Elastic::PAR::Coup2Pow* params) : params_(params) {}

void Mat::Elastic::Coup2Pow::AddStrainEnergy(double& psi, const Core::LinAlg::Matrix<3, 1>& prinv,
    const Core::LinAlg::Matrix<3, 1>& modinv, const Core::LinAlg::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const int d = params_->d_;

  // strain energy: Psi = C (II_{\boldsymbol{C}}-3)^D
  // add to overall strain energy
  psi += c * pow((prinv(1) - 3.), d);
}

void Mat::Elastic::Coup2Pow::add_derivatives_principal(Core::LinAlg::Matrix<3, 1>& dPI,
    Core::LinAlg::Matrix<6, 1>& ddPII, const Core::LinAlg::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  const double c = params_->c_;
  const int d = params_->d_;

  // If d<2 the material model is not stress free in the reference configuration
  if (d < 2)
    FOUR_C_THROW(
        "The Elast_Coup2Pow - material only works for positive integer exponents, which are larger "
        "than two.");

  dPI(1) += c * d * pow((prinv(1) - 3.), d - 1.);

  if (d == 2)
    ddPII(1) += (c * d * d - c * d);
  else
    ddPII(1) += (c * d * d - c * d) * pow((prinv(1) - 3.), d - 2.);
}
FOUR_C_NAMESPACE_CLOSE
