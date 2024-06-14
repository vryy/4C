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
      explicit IsotropicGrowthStrategy(const Core::Mat::PAR::Parameter::Data& matdata);

      std::unique_ptr<MIXTURE::MixtureGrowthStrategy> create_growth_strategy() override;
    };
  }  // namespace PAR

  /*!
   * @brief Growth is modeled as an inelastic volumentric expansion of the whole cell (isotropic).
   */
  class IsotropicGrowthStrategy : public MIXTURE::MixtureGrowthStrategy
  {
   public:
    [[nodiscard]] bool has_inelastic_growth_deformation_gradient() const override { return true; };

    void evaluate_inverse_growth_deformation_gradient(Core::LinAlg::Matrix<3, 3>& iFgM,
        const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
        int gp) const override;

    void evaluate_growth_stress_cmat(const MIXTURE::MixtureRule& mixtureRule,
        double currentReferenceGrowthScalar,
        const Core::LinAlg::Matrix<1, 6>& dCurrentReferenceGrowthScalarDC,
        const Core::LinAlg::Matrix<3, 3>& F, const Core::LinAlg::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>& S_stress,
        Core::LinAlg::Matrix<6, 6>& cmat, const int gp, const int eleGID) const override;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif