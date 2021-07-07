/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a penalty term like growth strategy

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef STIFFNESS_GROWTH_STRATEGY_H
#define STIFFNESS_GROWTH_STRATEGY_H

#include "mixture_growth_strategy.H"

namespace MIXTURE
{
  class StiffnessGrowthStrategy;

  namespace PAR
  {
    class StiffnessGrowthStrategy : public MIXTURE::PAR::MixtureGrowthStrategy
    {
     public:
      explicit StiffnessGrowthStrategy(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      std::unique_ptr<MIXTURE::MixtureGrowthStrategy> CreateGrowthStrategy() override;

      const double kappa_;
    };
  }  // namespace PAR

  /*!
   * @brief Growth modeled as an elastic expansion of the whole cell with a penalty type formulation
   *
   * This growth strategy uses a penalty term that ensures constant spatial density. The penalty
   * term has the following form:
   *
   * \f[
   *     \Psi = \frac{1}{2} \kappa (|F| - \frac{\rho_0(s)}{\rho_0(s=0)})^2
   * \f]
   *
   * The model is based on Braeu et al (2019),
   * https://link.springer.com/article/10.1007%2Fs10237-018-1084-x
   */
  class StiffnessGrowthStrategy : public MIXTURE::MixtureGrowthStrategy
  {
   public:
    explicit StiffnessGrowthStrategy(MIXTURE::PAR::StiffnessGrowthStrategy* params);

    void EvaluateInverseGrowthDeformationGradient(LINALG::Matrix<3, 3>& iFgM,
        const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
        int gp) const override;

    void AddGrowthStressCmat(const MIXTURE::MixtureRule& mixtureRule,
        double currentReferenceGrowthScalar, const LINALG::Matrix<3, 3>& F,
        const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
        LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, const int gp,
        const int eleGID) const override;

   private:
    ///! growth parameters as defined in the input file
    const PAR::StiffnessGrowthStrategy* params_{};
  };
}  // namespace MIXTURE

#endif  // STIFFNESS_GROWTH_STRATEGY_H