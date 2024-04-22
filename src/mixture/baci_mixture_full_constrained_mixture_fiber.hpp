/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a 1D remodel fiber
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_FULL_CONSTRAINED_MIXTURE_FIBER_HPP
#define FOUR_C_MIXTURE_FULL_CONSTRAINED_MIXTURE_FIBER_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_mixture_full_constrained_mixture_fiber_adaptive_history.hpp"
#include "baci_mixture_growth_evolution_linear_cauchy_poisson_turnover.hpp"

#include <Sacado.hpp>

#include <functional>
#include <limits>
#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::COMM
{
  class PackBuffer;
}
namespace MIXTURE
{
  template <typename T>
  class RemodelFiberMaterial;

  template <typename Number>
  struct MassIncrement
  {
    Number reference_stretch = 1.0;
    Number growth_scalar = 1.0;
    Number growth_scalar_production_rate = 0.0;
    double deposition_time = 0.0;
  };

  template <typename Number>
  bool IsAlmostEqual(
      const MassIncrement<Number>& inc1, const MassIncrement<Number>& inc2, const double tolerance)
  {
    return std::abs(inc1.reference_stretch - inc2.reference_stretch) < tolerance &&
           std::abs(inc1.growth_scalar - inc2.growth_scalar) < tolerance &&
           std::abs(inc1.growth_scalar_production_rate - inc2.growth_scalar_production_rate) <
               tolerance &&
           std::abs(inc1.deposition_time - inc2.deposition_time) < tolerance;
  }

  enum class HistoryAdaptionStrategy
  {
    none,
    window,
    model_equation,
    higher_order_integration
  };


  template <typename Number>
  struct State
  {
    Number lambda_f;
  };

  template <typename Number>
  struct DepositionHistoryInterval
  {
    std::vector<MassIncrement<Number>> timesteps = {};
    TimestepAdaptivityInfo adaptivity_info = {};
    double base_dt = 0.0;
  };

  template <typename Number>
  using DepositionHistory = std::vector<DepositionHistoryInterval<Number>>;

  /*!
   * @brief A full constrained mixture fiber based on the theory of Humphrey and Rajagopal (2002)
   * (https://doi.org/10.1142/S0218202502001714)
   *
   * This model generally assumes the deposition of new mass at every point with a specific
   * prestretch. Extant material is degraded over time with a Poisson degradation process. The
   * material itself stores the needed history variables. If activated, the number of history
   * variables is dynamically adapted to ensure efficient memory usage and fast evaluation times
   * while keeping the integration error low.
   *
   * @note This model is expensive in memory usage compared to the homogenized constrained mixture
   * fiber.
   *
   * @tparam Number The type of the number to be used. The default should be double, but also FAD
   * types can be used to check the exactness of the derivatives.
   */
  template <typename Number = double>
  class FullConstrainedMixtureFiber
  {
   public:
    FullConstrainedMixtureFiber(std::shared_ptr<const RemodelFiberMaterial<Number>> material,
        LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<Number> growth_evolution,
        Number lambda_pre, HistoryAdaptionStrategy adaptive_history_strategy, bool growth_enabled);

    /*!
     * @brief Pack all internal data into tha #data
     *
     * @param data (out) : buffer to serialize data to.
     */
    void Pack(CORE::COMM::PackBuffer& data) const;

    /*!
     * @brief Unpack all internal data that was previously packed by #Pack(CORE::COMM::PackBuffer&)
     *
     * @param position (in/out) : Position, where to start reading
     * @param data (in) : Vector of chars to extract data from
     */
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

    /// @brief Updates previous history data
    void Update();

    /*!
     * @brief Change the deposition stretch during simulation. This also changes the homeostatic
     * stretch of the fiber.
     *
     * @param lambda_pre (in) : Scalar for the prestretch of the fiber (in fiber direction)
     */
    void SetDepositionStretch(double lambda_pre);

    /*!
     * @brief Set deformation state of the fiber and recomputes all needed quantities for
     * evaluation
     *
     * @note This method has to be called before any Evaluation
     *
     * @param lambda_f (in) : total stretch in fiber direction
     * @param time (in) : Total time
     * @param dt (in) : Timestep
     */
    void RecomputeState(Number lambda_f, double time, double dt);

    /*!
     * @brief Reinitialize the history of the full constrained mixture fiber
     *
     * Call this method if a discontinuity occurs within lambda_f. The time must be the same as the
     * last inserted snapshot.
     *
     * @note This method must be called at least once at the beginning of a growth and remodeling
     * period.
     *
     * @param lambda_f (in) : Stretch of the fiber
     * @param time (in) : Total tiem
     */
    void ReinitializeHistory(Number lambda_f, double time);

    /*!
     * @brief Returns the last time record in the history, or 0.0 if it is empty
     *
     * @return double Last stored time in history
     */
    [[nodiscard]] double GetLastTimeInHistory() const;

    /*!
     * @brief Adds the time delta to all items in the history.
     *
     * @note This method is useful for a time interval without growth and remodeling.
     *
     * @param delta_time
     */
    void AddTime(double delta_time);

    ///@brief Integrands that need to be integrated over for the full constrained mixture model
    ///@{
    [[nodiscard]] Number GrowthScalarIntegrand(
        const MassIncrement<Number>& mass_increment, double time) const;
    [[nodiscard]] Number DGrowthScalarIntegrandDProductionRate(
        const MassIncrement<Number>& mass_increment, double time) const;
    [[nodiscard]] Number DGrowthScalarIntegrandDGrowthScalar(
        const MassIncrement<Number>& mass_increment, double time) const;
    [[nodiscard]] Number ScaledCauchyStressIntegrand(
        const MassIncrement<Number>& mass_increment, double time, Number current_lambda_f) const;
    [[nodiscard]] Number DScaledCauchyStressIntegrandDProductionRate(
        const MassIncrement<Number>& mass_increment, double time, Number current_lambda_f) const;
    [[nodiscard]] Number DScaledCauchyStressIntegrandDGrowthScalar(
        const MassIncrement<Number>& mass_increment, double time, Number current_lambda_f) const;
    [[nodiscard]] Number DScaledCauchyStressIntegrandDLambdaFSq(
        const MassIncrement<Number>& mass_increment, double time, Number current_lambda_f) const;
    [[nodiscard]] Number DScaledCauchyStressIntegrandDLambdaRefSq(
        const MassIncrement<Number>& mass_increment, double time, Number current_lambda_f) const;
    ///@}

    ///@brief Evaluator methods to compute the linearization of the local Newton method
    ///@{
    [[nodiscard]] Number EvaluateDResiduumGrowthScalarDLambdaFSq() const;
    [[nodiscard]] Number EvaluateDResiduumCauchyStressDLambdaFSq() const;
    ///@}

    /*!
     * @brief Returns an evaluator object that is passed to the local newton solver.
     *
     * The evaluator object returns an object that is callable with a CORE::LINALG::Matrix and
     * returns the residuum and the derivative.
     *
     * @return std::function<std::tuple<CORE::LINALG::Matrix<2, 1, Number>,
     * CORE::LINALG::Matrix<2, 2, Number>>(const CORE::LINALG::Matrix<2, 1, Number>&)>
     */
    [[nodiscard]] std::function<std::tuple<CORE::LINALG::Matrix<2, 1, Number>,
        CORE::LINALG::Matrix<2, 2, Number>>(const CORE::LINALG::Matrix<2, 1, Number>&)>
    GetLocalNewtonEvaluator() const;

    /*!
     * @brief Evaluate the Cauchy stress of the fiber by integration over the deposition history
     *
     * @param lambda_f Current fiber stretch
     * @return Number
     */
    [[nodiscard]] Number ComputeHistoryCauchyStress(Number lambda_f) const;

    [[nodiscard]] Number EvaluateLambdaRef(Number lambda_f) const;
    [[nodiscard]] Number EvaluateDLambdaRefSqDLambdaFSq(Number lambda_f) const;

    MassIncrement<Number> EvaluateCurrentMassIncrement(
        Number growth_scalar, Number cauchy_stress) const
    {
      const Number growth_scalar_production_rate =
          growth_evolution_.EvaluateTrueMassProductionRate((cauchy_stress - sig_h_) / sig_h_);
      return {EvaluateLambdaRef(current_state_.lambda_f), growth_scalar,
          growth_scalar_production_rate, current_time_};
    }

    void ComputeInternalVariables();

    /*!
     * @brief Return the current Cauchy fiber stress
     *
     * @return Number
     */
    [[nodiscard]] Number EvaluateCurrentCauchyStress() const { return computed_sigma_; }


    /*!
     * @brief Return the current 2. Piola Kirchhoff fiber stress
     *
     * @return Number
     */
    [[nodiscard]] Number EvaluateCurrentSecondPKStress() const
    {
      return computed_sigma_ / std::pow(current_state_.lambda_f, 2);
    }

    /*!
     * @brief Returns the derivative of the 2. Piola Kirchhoff fiber stress tensor w.r.t. the
     * squared fiber stretch
     *
     * @return Number
     */
    [[nodiscard]] Number EvaluateDCurrentFiberPK2StressDLambdafsq() const;

    /// homeostatic quantities
    /// @{
    Number sig_h_ = 0.0;
    Number lambda_pre_ = 1.0;
    /// @}

    /// states
    /// @{
    State<Number> current_state_ = {};
    /// @}

    /// Strain energy function of the fiber
    const std::shared_ptr<const RemodelFiberMaterial<Number>> fiber_material_;

    /// Growth evolution equation
    LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<Number> growth_evolution_;

    /// Flag whether growth is enabled
    bool growth_enabled_ = true;

    /// The deposition time of the initially present mass
    double reference_time_ = 0.0;

    /// A current time shit that is resetted after each timestep
    double current_time_shift_ = 0.0;

    HistoryAdaptionStrategy adaptive_history_strategy_ = HistoryAdaptionStrategy::none;
    DepositionHistory<Number> history_{};
    std::size_t window_size = 0;

    /// current data
    double current_time_ = 0.0;

    Number computed_growth_scalar_ = 1.0;
    Number computed_sigma_ = 0.0;
    Number computed_dgrowth_scalar_dlambda_f_sq_ = 0.0;
    Number computed_dsigma_dlambda_f_sq_ = 0.0;

    Number adaptive_tolerance_ = 1e-6;

#ifdef FOUR_C_ENABLE_ASSERTIONS
    bool state_is_set_ = false;
#endif
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif