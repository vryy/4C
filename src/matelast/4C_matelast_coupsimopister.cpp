/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the Simo and Pister material model except the volumetric term
\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupsimopister.hpp"

#include "4C_material_parameter_base.hpp"

#include <limits>

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::CoupSimoPister::CoupSimoPister(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata), mue_(matdata.parameters.Get<double>("MUE"))
{
}

Mat::Elastic::CoupSimoPister::CoupSimoPister(Mat::Elastic::PAR::CoupSimoPister* params)
    : params_(params)
{
}

void Mat::Elastic::CoupSimoPister::AddStrainEnergy(double& psi,
    const Core::LinAlg::Matrix<3, 1>& prinv, const Core::LinAlg::Matrix<3, 1>& modinv,
    const Core::LinAlg::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  // material Constant mu
  const double mue = params_->mue_;

  // strain energy: \Psi = 0.5*\mu(I_1-3) - \mu log(J)
  // add to overall strain energy
  psi += 0.5 * mue * (prinv(0) - 3.) - mue * log(std::pow(prinv(2), 0.5));
}

void Mat::Elastic::CoupSimoPister::add_derivatives_principal(Core::LinAlg::Matrix<3, 1>& dPI,
    Core::LinAlg::Matrix<6, 1>& ddPII, const Core::LinAlg::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  const double mue = params_->mue_;

  dPI(0) += 0.5 * mue;
  dPI(2) -= 0.5 * mue / prinv(2);

  ddPII(2) += 0.5 * mue / (prinv(2) * prinv(2));
}
FOUR_C_NAMESPACE_CLOSE
