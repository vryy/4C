/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an isotropic growth strategy for the growth remodel mixture rule

\level 3
*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_growth_strategy_isotropic.hpp"

#include "4C_mixture_growth_strategy.hpp"

FOUR_C_NAMESPACE_OPEN

MIXTURE::PAR::IsotropicGrowthStrategy::IsotropicGrowthStrategy(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : MIXTURE::PAR::MixtureGrowthStrategy(matdata)
{
}

std::unique_ptr<MIXTURE::MixtureGrowthStrategy>
MIXTURE::PAR::IsotropicGrowthStrategy::create_growth_strategy()
{
  return std::make_unique<MIXTURE::IsotropicGrowthStrategy>();
}

void MIXTURE::IsotropicGrowthStrategy::evaluate_inverse_growth_deformation_gradient(
    Core::LinAlg::Matrix<3, 3>& iFgM, const MIXTURE::MixtureRule& mixtureRule,
    double currentReferenceGrowthScalar, int gp) const
{
  iFgM.Clear();

  for (int i = 0; i < 3; ++i)
  {
    iFgM(i, i) = std::pow(currentReferenceGrowthScalar, -1.0 / 3.0);
  }
}

void MIXTURE::IsotropicGrowthStrategy::evaluate_growth_stress_cmat(
    const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
    const Core::LinAlg::Matrix<1, 6>& dCurrentReferenceGrowthScalarDC,
    const Core::LinAlg::Matrix<3, 3>& F, const Core::LinAlg::Matrix<6, 1>& E_strain,
    Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>& S_stress,
    Core::LinAlg::Matrix<6, 6>& cmat, const int gp, const int eleGID) const
{
  S_stress.Clear();
  cmat.Clear();
}
FOUR_C_NAMESPACE_CLOSE
