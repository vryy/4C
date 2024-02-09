/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a growth strategy for anisotropic growth

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_MIXTURE_GROWTH_STRATEGY_ANISOTROPIC_HPP
#define BACI_MIXTURE_GROWTH_STRATEGY_ANISOTROPIC_HPP

#include "baci_config.hpp"

#include "baci_mat_anisotropy_extension_default.hpp"
#include "baci_mixture_growth_strategy.hpp"

BACI_NAMESPACE_OPEN

namespace MIXTURE
{
  class AnisotropicGrowthStrategy;

  namespace PAR
  {
    class AnisotropicGrowthStrategy : public MIXTURE::PAR::MixtureGrowthStrategy
    {
     public:
      explicit AnisotropicGrowthStrategy(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      std::unique_ptr<MIXTURE::MixtureGrowthStrategy> CreateGrowthStrategy() override;

      const int init_mode_;
      const int fiber_id_;

      /// structural tensor strategy
      Teuchos::RCP<MAT::ELASTIC::StructuralTensorStrategyBase> structural_tensor_strategy_;
    };
  }  // namespace PAR

  /*!
   * @brief Growth is modeled as an inelastic expansion of the whole cell in one predefined
   * direction
   *
   * The direction of growth can be specified with a fiber with the fiber id specified in FIBER_ID
   *
   */
  class AnisotropicGrowthStrategy : public MIXTURE::MixtureGrowthStrategy
  {
   public:
    explicit AnisotropicGrowthStrategy(MIXTURE::PAR::AnisotropicGrowthStrategy* params);

    void PackMixtureGrowthStrategy(CORE::COMM::PackBuffer& data) const override;

    void UnpackMixtureGrowthStrategy(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy) override;

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

   private:
    ///! growth parameters as defined in the input file
    const PAR::AnisotropicGrowthStrategy* params_{};

    /// Anisotropy extension that manages fibers and structural tensors
    MAT::DefaultAnisotropyExtension<1> anisotropyExtension_;
  };
}  // namespace MIXTURE

BACI_NAMESPACE_CLOSE

#endif  // MIXTURE_GROWTH_STRATEGY_ANISOTROPIC_H