/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an isotropic general power-type material in terms of the first Cauchy-Green
invariant

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coup1pow.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::Elastic::PAR::Coup1Pow::Coup1Pow(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      c_(matdata.parameters.Get<double>("C")),
      d_(matdata.parameters.Get<int>("D"))
{
}

Mat::Elastic::Coup1Pow::Coup1Pow(Mat::Elastic::PAR::Coup1Pow* params) : params_(params) {}

void Mat::Elastic::Coup1Pow::AddStrainEnergy(double& psi, const Core::LinAlg::Matrix<3, 1>& prinv,
    const Core::LinAlg::Matrix<3, 1>& modinv, const Core::LinAlg::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const int d = params_->d_;

  // strain energy: Psi = C (I_{\boldsymbol{C}}-3)^D
  // add to overall strain energy
  psi += c * pow((prinv(0) - 3.), d);
}

void Mat::Elastic::Coup1Pow::add_derivatives_principal(Core::LinAlg::Matrix<3, 1>& dPI,
    Core::LinAlg::Matrix<6, 1>& ddPII, const Core::LinAlg::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  const double c = params_->c_;
  const double d = params_->d_;

  // If d<2 the material model is not stress free in the reference configuration
  if (d < 2)
    FOUR_C_THROW(
        "The Elast_Coup1Pow - material only works for positive integer exponents, which are larger "
        "than two.");

  dPI(0) += c * d * pow((prinv(0) - 3.), d - 1.);

  if (d == 2)
    ddPII(0) += (c * d * d - c * d);
  else
    ddPII(0) += (c * d * d - c * d) * pow((prinv(0) - 3.), d - 2.);
}
FOUR_C_NAMESPACE_CLOSE
