/*----------------------------------------------------------------------*/
/*! \file

\brief 1D constitutive equation for growth (Linear in Cauchy fiber stress)

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_GROWTH_EVOLUTION_LINEAR_CAUCHY_POISSON_TURNOVER_HPP
#define FOUR_C_MIXTURE_GROWTH_EVOLUTION_LINEAR_CAUCHY_POISSON_TURNOVER_HPP


#include "baci_config.hpp"

#include <cmath>

FOUR_C_NAMESPACE_OPEN

namespace MIXTURE
{
  /*!
   * @brief Defines the evolution equation of net 1D mass production during growth and remodeling
   */
  template <typename Number>
  class LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution
  {
   public:
    /*!
     * \brief Default constructor
     *
     * \param k_sig Gain parameter of linear growth equation
     */
    LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution(double k_sig, double decay_time)
        : k_sig_(k_sig), decay_time_(decay_time)
    {
    }

    [[nodiscard]] Number EvaluateTrueMassProductionRate(Number delta_sig) const
    {
      return k_sig_ * delta_sig + 1.0 / decay_time_;
    }

    [[nodiscard]] Number EvaluateTrueMassRemovalRate(Number delta_sig) const
    {
      return -1.0 / decay_time_;
    }

    [[nodiscard]] Number EvaluateDTrueMassProductionRateDSig(Number delta_sig) const
    {
      return k_sig_;
    }

    [[nodiscard]] Number EvaluateDTrueMassRemovalRateDSig(Number delta_sig) const { return 0.0; }

    [[nodiscard]] inline double EvaluateSurvivalFunction(double delta_time) const
    {
      return std::exp(-delta_time / decay_time_);
    }

    [[nodiscard]] inline double EvaluateSurvivalFunctionIntegration(
        const double current_time, const double t1, const double t2) const
    {
      return decay_time_ * (std::exp(-(current_time - t2) / decay_time_) -
                               std::exp(-(current_time - t1) / decay_time_));
    }

    /// Gain-Parameter controlling the speed of the growth process
    double k_sig_;

    /// decay time controlling the speed of turnover
    double decay_time_;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif
