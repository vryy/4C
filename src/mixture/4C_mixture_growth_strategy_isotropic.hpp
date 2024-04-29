/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of an isotropic growth strategy for the growth remodel mixture rule

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MIXTURE_GROWTH_STRATEGY_ISOTROPIC_HPP
#define FOUR_C_MIXTURE_GROWTH_STRATEGY_ISOTROPIC_HPP

#include "4C_config.hpp"

#include "4C_mixture_growth_strategy.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MIXTURE
{
  class IsotropicGrowthStrategy;

  namespace PAR
  {
    class IsotropicGrowthStrategy : public MIXTURE::PAR::MixtureGrowthStrategy
    {
     public:
      explicit IsotropicGrowthStrategy(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

      std::unique_ptr<MIXTURE::MixtureGrowthStrategy> CreateGrowthStrategy() override;
    };
  }  // namespace PAR

  /*!
   * @brief Growth is modeled as an inelastic volumentric expansion of the whole cell (isotropic).
   */
  class IsotropicGrowthStrategy : public MIXTURE::MixtureGrowthStrategy
  {
   public:
    [[nodiscard]] bool HasInelasticGrowthDeformationGradient() const override { return true; };

    void EvaluateInverseGrowthDeformationGradient(CORE::LINALG::Matrix<3, 3>& iFgM,
        const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
        int gp) const override;

    void EvaluateGrowthStressCmat(const MIXTURE::MixtureRule& mixtureRule,
        double currentReferenceGrowthScalar,
        const CORE::LINALG::Matrix<1, 6>& dCurrentReferenceGrowthScalarDC,
        const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, const int gp, const int eleGID) const override;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif