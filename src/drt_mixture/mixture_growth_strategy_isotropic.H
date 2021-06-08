/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of an isotropic growth strategy for the growth remodel mixture rule

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef ISOTROPIC_GROWTH_STRATEGY_H
#define ISOTROPIC_GROWTH_STRATEGY_H

#include "mixture_growth_strategy.H"

namespace MIXTURE
{
  class IsotropicGrowthStrategy;

  namespace PAR
  {
    class IsotropicGrowthStrategy : public MIXTURE::PAR::MixtureGrowthStrategy
    {
     public:
      explicit IsotropicGrowthStrategy(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      std::unique_ptr<MIXTURE::MixtureGrowthStrategy> CreateGrowthStrategy() override;
    };
  }  // namespace PAR

  class IsotropicGrowthStrategy : public MIXTURE::MixtureGrowthStrategy
  {
   public:
    void EvaluateInverseGrowthDeformationGradient(LINALG::Matrix<3, 3>& iFgM,
        const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
        int gp) const override;

    void AddGrowthStressCmat(const MIXTURE::MixtureRule& mixtureRule,
        double currentReferenceGrowthScalar, const LINALG::Matrix<3, 3>& F,
        const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
        LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, const int gp,
        const int eleGID) const override;
  };
}  // namespace MIXTURE

#endif  // ISOTROPIC_GROWTH_STRATEGY_H