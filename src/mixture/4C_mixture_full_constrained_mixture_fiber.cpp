/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a 1D full constrained mixture fiber
\level 3
*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_full_constrained_mixture_fiber.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_mixture_constituent_remodelfiber_material.hpp"
#include "4C_mixture_full_constrained_mixture_fiber_adaptive_history.hpp"
#include "4C_mixture_growth_evolution_linear_cauchy_poisson_turnover.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_local_integration.hpp"
#include "4C_utils_local_newton.hpp"

#include <Sacado.hpp>

#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <type_traits>

FOUR_C_NAMESPACE_OPEN

namespace
{
  bool IsNear(const double value1, const double value2, const double epsilon = 1e-15)
  {
    if (std::abs(value1 - value2) <= epsilon) return true;
    return std::abs(value1 - value2) <= epsilon * std::max(std::abs(value1), std::abs(value2));
  }

  template <typename Number>
  static inline Number ComputeAbsoluteErrorWhenSkippingSnapshot(
      std::size_t i, const std::vector<std::tuple<double, Number>>& evaluated_integrand)
  {
    FOUR_C_ASSERT(i > 0,
        "This function should only be called for an item in the interior such that it can still be "
        "integrated with a Simpson's rule if we remove point i (0 < i < size-2)");
    FOUR_C_ASSERT(i < evaluated_integrand.size() - 2,
        "This function should only be called for an item in the interior such that it can still be "
        "integrated with a Simpson's rule if we remove point i (0 < i < size-2)");

    Number full_integration =
        Core::UTILS::IntegrateSimpsonStep<std::tuple<double, Number>>(
            evaluated_integrand[i - 1], evaluated_integrand[i], evaluated_integrand[i + 1]) +
        Core::UTILS::IntegrateSimpsonStepBC<std::tuple<double, Number>>(
            evaluated_integrand[i + 0], evaluated_integrand[i + 1], evaluated_integrand[i + 2]);
    Number skipped_integration = Core::UTILS::IntegrateSimpsonStep<std::tuple<double, Number>>(
        evaluated_integrand[i - 1], evaluated_integrand[i + 1], evaluated_integrand[i + 2]);

    return std::abs(full_integration - skipped_integration);
  }

  template <typename Number>
  static inline std::vector<bool> GetEraseItemsBasedOnRelativeError(
      const std::vector<Number>& rel_error, Number tolerance)
  {
    std::vector<std::size_t> unit(rel_error.size() - 3);
    std::iota(unit.begin(), unit.end(), 1);

    auto num_items_below_tolerance = unit.begin();
    std::for_each(rel_error.begin(), rel_error.end(),
        [&](Number item) { num_items_below_tolerance += item < tolerance; });

    std::partial_sort(unit.begin(), num_items_below_tolerance, unit.end(),
        [&](std::size_t a, std::size_t b) { return rel_error[a] < rel_error[b]; });

    std::vector<bool> erase_vector(rel_error.size(), false);

    Number total_relative_error = 0.0;
    for (std::size_t index : unit)
    {
      // ensure no adjacent deletions
      if (erase_vector[index - 1] || erase_vector[index + 1]) continue;
      total_relative_error += rel_error[index];
      if (total_relative_error < tolerance) erase_vector[index] = true;
    }

    return erase_vector;
  }

  template <typename Number, typename Integrand>
  static inline Number IntegrateOverDepositionHistory(
      const MIXTURE::DepositionHistory<Number>& history, Integrand integrand)
  {
    Number integration_result = 0;
    for (const auto& interval : history)
    {
      integration_result += Core::UTILS::IntegrateSimpsonTrapezoidal(interval.timesteps,
          [&](const MIXTURE::MassIncrement<Number>& increment)
          { return std::make_tuple(increment.deposition_time, integrand(increment)); });
    }

    return integration_result;
  }

  template <typename Number, typename Integrand>
  static inline std::tuple<Number, Number> IntegrateLastTimestepWithDerivative(
      const MIXTURE::DepositionHistoryInterval<Number>& interval,
      const MIXTURE::MassIncrement<Number>& current_increment, Integrand integrand,
      Number history_integration)
  {
    const size_t size = interval.timesteps.size();
    if (size == 0) FOUR_C_THROW("The history is empty. I cannot integrate.");
    if (size == 1)
    {
      // can only apply trapezoidal rule
      const auto [integration, derivative] =
          Core::UTILS::IntegrateTrapezoidalStepAndReturnDerivativeB<std::tuple<double, Number>>(
              {interval.timesteps[0].deposition_time, integrand(interval.timesteps[0])},
              {current_increment.deposition_time, integrand(current_increment)});
      return std::make_tuple(integration + history_integration, derivative);
    }

    const auto [integration, derivative] =
        Core::UTILS::IntegrateSimpsonStepBCAndReturnDerivativeC<std::tuple<double, Number>>(
            {interval.timesteps[size - 2].deposition_time, integrand(interval.timesteps[size - 2])},
            {interval.timesteps[size - 1].deposition_time, integrand(interval.timesteps[size - 1])},
            {current_increment.deposition_time, integrand(current_increment)});

    return std::make_tuple(integration + history_integration, derivative);
  }

  template <typename Number>
  static inline Number EvaluateI4(Number lambda_e)
  {
    return std::pow(lambda_e, 2);
  }

  template <typename Number>
  static inline Number EvaluateDI4DlLambdaESq(Number lambda_e)
  {
    return 1.0;
  }

  template <typename Number>
  static inline Number EvaluateFiberMaterialCauchyStress(
      const MIXTURE::RemodelFiberMaterial<Number>& fiber_material, Number lambda_e)
  {
    const Number I4 = EvaluateI4(lambda_e);
    return fiber_material.GetCauchyStress(I4);
  }

  template <typename Number>
  static inline Number EvaluateDFiberMaterialCauchyStressDLambdaESq(
      const MIXTURE::RemodelFiberMaterial<Number>& fiber_material, Number lambda_e)
  {
    const Number I4 = EvaluateI4(lambda_e);
    const Number DI4DLambdaESq = EvaluateDI4DlLambdaESq(lambda_e);
    return fiber_material.GetDCauchyStressDI4(I4) * DI4DLambdaESq;
  }

  template <typename Number>
  static inline void ReinitializeState(
      MIXTURE::FullConstrainedMixtureFiber<Number>& fiber, const Number lambda_f, const double time)
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    fiber.state_is_set_ = true;
#endif
    fiber.current_state_.lambda_f = lambda_f;
    fiber.current_time_ = time;
  }

  template <typename Number>
  void UpdateBaseDeltaTime(
      MIXTURE::DepositionHistoryInterval<Number>& deposition_history_inverval, const double dt)
  {
    if (deposition_history_inverval.base_dt <= 0)
    {
      // timestep is set for the first time
      deposition_history_inverval.base_dt = dt;
      return;
    }

    if (!IsNear(deposition_history_inverval.base_dt, dt))
    {
      FOUR_C_THROW(
          "The timestep is not constant within the interval. The interval currently relies on a "
          "constant timestep of %f. You are stepping with %f (err = %f). You need to extend the "
          "implementation such that it can also handle adaptive/non equidistant timestepping.",
          deposition_history_inverval.base_dt, dt,
          std::abs(deposition_history_inverval.base_dt - dt));
    }
  }

  template <typename Number>
  Number IntegrateBooleStep(const std::array<std::tuple<double, Number>, 5>& data)
  {
    const double interval_size = std::get<0>(data[4]) - std::get<0>(data[0]);

    return 1.0 / 90.0 * interval_size *
           (7 * std::get<1>(data[0]) + 32 * std::get<1>(data[1]) + 12 * std::get<1>(data[2]) +
               32 * std::get<1>(data[3]) + 7 * std::get<1>(data[4]));
  }
}  // namespace

template <typename Number>
MIXTURE::FullConstrainedMixtureFiber<Number>::FullConstrainedMixtureFiber(
    std::shared_ptr<const MIXTURE::RemodelFiberMaterial<Number>> material,
    MIXTURE::LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<Number> growth_evolution,
    Number lambda_pre, HistoryAdaptionStrategy adaptive_history_strategy, bool growth_enabled)
    : lambda_pre_(lambda_pre),
      fiber_material_(std::move(material)),
      growth_evolution_(growth_evolution),
      growth_enabled_(growth_enabled),
      adaptive_history_strategy_(adaptive_history_strategy),
      current_time_(0),
      computed_growth_scalar_(1.0)
{
  sig_h_ = EvaluateFiberMaterialCauchyStress<Number>(*fiber_material_, lambda_pre_);
  computed_sigma_ = sig_h_;
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::Pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("Packing and Unpacking is currently only implemented for the double-specialization");
}

template <>
void MIXTURE::FullConstrainedMixtureFiber<double>::Pack(Core::Communication::PackBuffer& data) const
{
  data.AddtoPack(sig_h_);
  data.AddtoPack(lambda_pre_);

  data.AddtoPack(current_state_.lambda_f);

  data.AddtoPack(reference_time_);
  data.AddtoPack(current_time_shift_);

  data.AddtoPack(history_.size());
  for (const auto& interval : history_)
  {
    data.AddtoPack(interval.timesteps.size());
    for (const auto& item : interval.timesteps)
    {
      data.AddtoPack(item.reference_stretch);
      data.AddtoPack(item.growth_scalar);
      data.AddtoPack(item.growth_scalar_production_rate);
      data.AddtoPack(item.deposition_time);
    }

    data.AddtoPack(interval.base_dt);
    interval.adaptivity_info.Pack(data);
  }

  data.AddtoPack(current_time_);

  data.AddtoPack(computed_growth_scalar_);
  data.AddtoPack(computed_sigma_);
  data.AddtoPack(computed_dgrowth_scalar_dlambda_f_sq_);
  data.AddtoPack(computed_dsigma_dlambda_f_sq_);
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  FOUR_C_THROW("Packing and Unpacking is currently only implemented for the double-specialization");
}

template <>
void MIXTURE::FullConstrainedMixtureFiber<double>::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  Core::Communication::ParObject::ExtractfromPack(position, data, sig_h_);
  Core::Communication::ParObject::ExtractfromPack(position, data, lambda_pre_);
  Core::Communication::ParObject::ExtractfromPack(position, data, current_state_.lambda_f);

  Core::Communication::ParObject::ExtractfromPack(position, data, reference_time_);
  Core::Communication::ParObject::ExtractfromPack(position, data, current_time_shift_);

  std::size_t size_of_history;
  Core::Communication::ParObject::ExtractfromPack(position, data, size_of_history);
  history_.resize(size_of_history);

  for (auto& interval : history_)
  {
    std::size_t size_of_interval;
    Core::Communication::ParObject::ExtractfromPack(position, data, size_of_interval);
    interval.timesteps.resize(size_of_interval);
    for (auto& item : interval.timesteps)
    {
      Core::Communication::ParObject::ExtractfromPack(position, data, item.reference_stretch);
      Core::Communication::ParObject::ExtractfromPack(position, data, item.growth_scalar);
      Core::Communication::ParObject::ExtractfromPack(
          position, data, item.growth_scalar_production_rate);
      Core::Communication::ParObject::ExtractfromPack(position, data, item.deposition_time);
    }

    Core::Communication::ParObject::ExtractfromPack(position, data, interval.base_dt);

    interval.adaptivity_info.Unpack(position, data);
  }


  Core::Communication::ParObject::ExtractfromPack(position, data, current_time_);


  Core::Communication::ParObject::ExtractfromPack(position, data, computed_growth_scalar_);
  Core::Communication::ParObject::ExtractfromPack(position, data, computed_sigma_);
  Core::Communication::ParObject::ExtractfromPack(
      position, data, computed_dgrowth_scalar_dlambda_f_sq_);
  Core::Communication::ParObject::ExtractfromPack(position, data, computed_dsigma_dlambda_f_sq_);
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::RecomputeState(
    const Number lambda_f, const double time, const double dt)
{
  ReinitializeState(*this, lambda_f, time);

  if (!growth_enabled_)
  {
    current_time_shift_ = get_last_time_in_history() - current_time_;
  }
  else
  {
    current_time_shift_ = 0.0;
  }

  if (growth_enabled_ && history_.size() > 0) UpdateBaseDeltaTime(history_.back(), dt);
  compute_internal_variables();
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::compute_history_cauchy_stress(
    const Number lambda_f) const
{
  const Number growth_scalar = std::invoke(
      [&]()
      {
        if (history_.size() > 0) return history_.back().timesteps.back().growth_scalar;
        return Number(1.0);
      });
  const double time = std::invoke(
      [&]()
      {
        if (history_.size() > 0) return history_.back().timesteps.back().deposition_time;
        return 0.0;
      });


  const auto current_scaled_cauchy_stress_integrand =
      [&](const MassIncrement<Number>& mass_increment)
  { return scaled_cauchy_stress_integrand(mass_increment, time, lambda_f); };

  const Number scaled_cauchy_stress_history =
      growth_evolution_.evaluate_survival_function(time) *
          EvaluateFiberMaterialCauchyStress<Number>(*fiber_material_, lambda_pre_ * lambda_f) +
      IntegrateOverDepositionHistory<Number>(history_, current_scaled_cauchy_stress_integrand);

  return scaled_cauchy_stress_history / growth_scalar;
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::growth_scalar_integrand(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time) const
{
  return mass_increment.growth_scalar_production_rate * mass_increment.growth_scalar *
         growth_evolution_.evaluate_survival_function(time - mass_increment.deposition_time);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::d_growth_scalar_integrand_d_production_rate(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time) const
{
  return mass_increment.growth_scalar *
         growth_evolution_.evaluate_survival_function(time - mass_increment.deposition_time);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::d_growth_scalar_integrand_d_growth_scalar(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time) const
{
  return mass_increment.growth_scalar_production_rate *
         growth_evolution_.evaluate_survival_function(time - mass_increment.deposition_time);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::scaled_cauchy_stress_integrand(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  return mass_increment.growth_scalar_production_rate * mass_increment.growth_scalar *
         growth_evolution_.evaluate_survival_function(time - mass_increment.deposition_time) *
         EvaluateFiberMaterialCauchyStress<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f);
}

template <typename Number>
Number
MIXTURE::FullConstrainedMixtureFiber<Number>::d_scaled_cauchy_stress_integrand_d_production_rate(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  FOUR_C_ASSERT(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  return mass_increment.growth_scalar *
         growth_evolution_.evaluate_survival_function(time - mass_increment.deposition_time) *
         EvaluateFiberMaterialCauchyStress<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f);
}

template <typename Number>
Number
MIXTURE::FullConstrainedMixtureFiber<Number>::d_scaled_cauchy_stress_integrand_d_growth_scalar(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  FOUR_C_ASSERT(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  return mass_increment.growth_scalar_production_rate *
         growth_evolution_.evaluate_survival_function(time - mass_increment.deposition_time) *
         EvaluateFiberMaterialCauchyStress<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::d_scaled_cauchy_stress_integrand_d_lambda_f_sq(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  FOUR_C_ASSERT(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  return mass_increment.growth_scalar_production_rate * mass_increment.growth_scalar *
         growth_evolution_.evaluate_survival_function(time - mass_increment.deposition_time) *
         EvaluateDFiberMaterialCauchyStressDLambdaESq<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f) *
         mass_increment.reference_stretch * mass_increment.reference_stretch;
}

template <typename Number>
Number
MIXTURE::FullConstrainedMixtureFiber<Number>::d_scaled_cauchy_stress_integrand_d_lambda_ref_sq(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  FOUR_C_ASSERT(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  return mass_increment.growth_scalar_production_rate * mass_increment.growth_scalar *
         growth_evolution_.evaluate_survival_function(time - mass_increment.deposition_time) *
         EvaluateDFiberMaterialCauchyStressDLambdaESq<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f) *
         current_lambda_f * current_lambda_f;
}

template <typename Number>
std::function<std::tuple<Core::LinAlg::Matrix<2, 1, Number>, Core::LinAlg::Matrix<2, 2, Number>>(
    const Core::LinAlg::Matrix<2, 1, Number>&)>
MIXTURE::FullConstrainedMixtureFiber<Number>::get_local_newton_evaluator() const
{
  FOUR_C_ASSERT(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  FOUR_C_ASSERT(current_time_shift_ == 0.0, "The timeshift should be zero if growth is enabled");
  FOUR_C_ASSERT(
      history_.size() > 0, "The history is empty. Please initialize it with ReinitializeHistory()");

  const auto current_growth_scalar_integrand = [&](const MassIncrement<Number>& mass_increment)
  { return growth_scalar_integrand(mass_increment, current_time_); };

  const auto current_dgrowth_scalar_integrand_dproduction_rate =
      [&](const MassIncrement<Number>& mass_increment)
  { return d_growth_scalar_integrand_d_production_rate(mass_increment, current_time_); };

  const auto current_dgrowth_scalar_integrand_dgrowth_scalar =
      [&](const MassIncrement<Number>& mass_increment)
  { return d_growth_scalar_integrand_d_growth_scalar(mass_increment, current_time_); };

  const auto current_scaled_cauchy_stress_integrand =
      [&](const MassIncrement<Number>& mass_increment)
  {
    return scaled_cauchy_stress_integrand(mass_increment, current_time_, current_state_.lambda_f);
  };

  const auto current_dscaled_cauchy_stress_integrand_dproduction_rate =
      [&](const MassIncrement<Number>& mass_increment)
  {
    return d_scaled_cauchy_stress_integrand_d_production_rate(
        mass_increment, current_time_, current_state_.lambda_f);
  };

  const auto current_dscaled_cauchy_stress_integrand_dgrowth_scalar =
      [&](const MassIncrement<Number>& mass_increment)
  {
    return d_scaled_cauchy_stress_integrand_d_growth_scalar(
        mass_increment, current_time_, current_state_.lambda_f);
  };

  // evaluate history quantities (except last timestep, which depends on the current Cauchy
  // stress)
  const Number growth_scalar_history =
      growth_evolution_.evaluate_survival_function(current_time_ - reference_time_) +
      IntegrateOverDepositionHistory<Number>(history_, current_growth_scalar_integrand);

  const Number scaled_cauchy_stress_history =
      growth_evolution_.evaluate_survival_function(current_time_ - reference_time_) *
          EvaluateFiberMaterialCauchyStress<Number>(
              *fiber_material_, lambda_pre_ * current_state_.lambda_f) +
      IntegrateOverDepositionHistory<Number>(history_, current_scaled_cauchy_stress_integrand);

  return [=](const Core::LinAlg::Matrix<2, 1, Number>& growth_scalar_and_cauchy_stress)
  {
    const Number growth_scalar = growth_scalar_and_cauchy_stress(0);
    const Number cauchy_stress = growth_scalar_and_cauchy_stress(1);

    const MassIncrement<Number> current_increment =
        evaluate_current_mass_increment(growth_scalar, cauchy_stress);

    const Number dmy_growth_scalar_production_rate_dsig =
        growth_evolution_.evaluate_d_true_mass_production_rate_d_sig(
            (cauchy_stress - sig_h_) / sig_h_) /
        sig_h_;
    const auto [my_growth_scalar, dmy_growth_scalar_dintegrand] =
        IntegrateLastTimestepWithDerivative<Number>(history_.back(), current_increment,
            current_growth_scalar_integrand, growth_scalar_history);

    const Number dmy_growth_scalar_dsig =
        dmy_growth_scalar_dintegrand *
        current_dgrowth_scalar_integrand_dproduction_rate(current_increment) *
        dmy_growth_scalar_production_rate_dsig;
    const Number dmy_growth_scalar_dgrowth_scalar =
        dmy_growth_scalar_dintegrand *
        current_dgrowth_scalar_integrand_dgrowth_scalar(current_increment);

    const auto [my_scaled_cauchy_stress, dmy_scaled_cauchy_stress_dintegrand] =
        IntegrateLastTimestepWithDerivative<Number>(history_.back(), current_increment,
            current_scaled_cauchy_stress_integrand, scaled_cauchy_stress_history);


    const Number dmy_scaled_cauchy_stress_dsig =
        dmy_scaled_cauchy_stress_dintegrand *
        current_dscaled_cauchy_stress_integrand_dproduction_rate(current_increment) *
        dmy_growth_scalar_production_rate_dsig;
    const Number dmy_scaled_cauchy_stress_dgrowth_scalar =
        dmy_scaled_cauchy_stress_dintegrand *
        current_dscaled_cauchy_stress_integrand_dgrowth_scalar(current_increment);

    Core::LinAlg::Matrix<2, 1, Number> residua;
    residua(0) = my_growth_scalar - growth_scalar;
    residua(1) = my_scaled_cauchy_stress / growth_scalar - cauchy_stress;

    Core::LinAlg::Matrix<2, 2, Number> derivative;
    derivative(0, 0) =
        dmy_growth_scalar_dgrowth_scalar - 1.0;  // d residuum growth scalar d growth scalar
    derivative(0, 1) = dmy_growth_scalar_dsig;   // d residuum growth scalar d sig
    derivative(1, 0) =
        (dmy_scaled_cauchy_stress_dgrowth_scalar * growth_scalar - my_scaled_cauchy_stress) /
        std::pow(growth_scalar, 2);  // d residuum sigma d growth scalar
    derivative(1, 1) =
        dmy_scaled_cauchy_stress_dsig / growth_scalar - 1.0;  // d residuum sigma d sigma

    return std::make_tuple(residua, derivative);
  };
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<
    Number>::evaluate_d_residuum_growth_scalar_d_lambda_f_sq() const
{
  return 0.0;
}


template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<
    Number>::evaluate_d_residuum_cauchy_stress_d_lambda_f_sq() const
{
  FOUR_C_ASSERT(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  FOUR_C_ASSERT(current_time_shift_ == 0.0, "The timeshift should be zero if growth is enabled");
  FOUR_C_ASSERT(
      history_.size() > 0, "The history is empty. Please initialize it with ReinitializeHistory()");

  const auto dscaled_cauchy_stress_integrand_d_lambda_f_sq =
      [&](const MassIncrement<Number>& mass_increment)
  {
    return d_scaled_cauchy_stress_integrand_d_lambda_f_sq(
        mass_increment, current_time_, current_state_.lambda_f);
  };
  const Number Dscaled_cauchy_stress_D_lambda_f_sq_history =
      growth_evolution_.evaluate_survival_function(current_time_ - reference_time_) *
          EvaluateDFiberMaterialCauchyStressDLambdaESq<Number>(
              *fiber_material_, lambda_pre_ * current_state_.lambda_f) *
          std::pow(lambda_pre_, 2) +
      IntegrateOverDepositionHistory<Number>(
          history_, dscaled_cauchy_stress_integrand_d_lambda_f_sq);

  const MassIncrement<Number> current_increment =
      evaluate_current_mass_increment(computed_growth_scalar_, computed_sigma_);

  const auto [dscaled_cauchy_stress_lambda_f_sq,
      ddscaled_cauchy_stress_integrand_d_lambda_f_sq_d_integrand] =
      IntegrateLastTimestepWithDerivative(history_.back(), current_increment,
          dscaled_cauchy_stress_integrand_d_lambda_f_sq,
          Dscaled_cauchy_stress_D_lambda_f_sq_history);


  return (dscaled_cauchy_stress_lambda_f_sq +
             ddscaled_cauchy_stress_integrand_d_lambda_f_sq_d_integrand *
                 d_scaled_cauchy_stress_integrand_d_lambda_ref_sq(
                     current_increment, current_time_, current_state_.lambda_f) *
                 evaluate_d_lambda_ref_sq_d_lambda_f_sq(current_state_.lambda_f)) /
         computed_growth_scalar_;
}


template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::EvaluateLambdaRef(Number lambda_f) const
{
  return lambda_pre_ / lambda_f;
}


template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::evaluate_d_lambda_ref_sq_d_lambda_f_sq(
    Number lambda_f) const
{
  return -std::pow(lambda_pre_, 2) / std::pow(lambda_f, 4);
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::compute_internal_variables()
{
  if (!growth_enabled_)
  {
    // in this case, I don't need to solve a local Newton method. I just need to do the integration
    // over the history until the moment.
    computed_dgrowth_scalar_dlambda_f_sq_ = 0;

    const auto current_growth_scalar_integrand = [&](const MassIncrement<Number>& mass_increment)
    { return growth_scalar_integrand(mass_increment, current_time_ + current_time_shift_); };

    const auto current_scaled_cauchy_stress_integrand =
        [&](const MassIncrement<Number>& mass_increment)
    {
      return scaled_cauchy_stress_integrand(
          mass_increment, current_time_ + current_time_shift_, current_state_.lambda_f);
    };

    const auto current_dscaled_cauchy_stress_integrand_d_lambda_f_sq =
        [&](const MassIncrement<Number>& mass_increment)
    {
      return d_scaled_cauchy_stress_integrand_d_lambda_f_sq(
          mass_increment, current_time_, current_state_.lambda_f);
    };

    computed_growth_scalar_ =
        growth_evolution_.evaluate_survival_function(
            current_time_ + current_time_shift_ - reference_time_) +
        IntegrateOverDepositionHistory<Number>(history_, current_growth_scalar_integrand);

    computed_sigma_ = (growth_evolution_.evaluate_survival_function(
                           current_time_ + current_time_shift_ - reference_time_) *
                              EvaluateFiberMaterialCauchyStress<Number>(
                                  *fiber_material_, lambda_pre_ * current_state_.lambda_f) +
                          IntegrateOverDepositionHistory<Number>(
                              history_, current_scaled_cauchy_stress_integrand)) /
                      computed_growth_scalar_;


    computed_dsigma_dlambda_f_sq_ =
        (growth_evolution_.evaluate_survival_function(
             current_time_ + current_time_shift_ - reference_time_) *
                EvaluateDFiberMaterialCauchyStressDLambdaESq<Number>(
                    *fiber_material_, lambda_pre_ * current_state_.lambda_f) *
                std::pow(lambda_pre_, 2) +
            IntegrateOverDepositionHistory<Number>(
                history_, current_dscaled_cauchy_stress_integrand_d_lambda_f_sq)) /
        computed_growth_scalar_;

    return;
  }

  //  Evaluate local Newton system to obtain the current Cauchy stress
  const auto EvaluateCurrentLocalNewtonLinearSystem = get_local_newton_evaluator();

  constexpr auto tolerance = 1e-10;
  constexpr auto max_iterations = 500;

  Core::LinAlg::Matrix<2, 1, Number> initial_guess;
  initial_guess(0) = computed_growth_scalar_;
  initial_guess(1) = computed_sigma_;
  auto [growth_scalar_and_sigma, K] = Core::UTILS::SolveLocalNewtonAndReturnJacobian(
      EvaluateCurrentLocalNewtonLinearSystem, initial_guess, tolerance, max_iterations);

  computed_growth_scalar_ = growth_scalar_and_sigma(0);
  computed_sigma_ = growth_scalar_and_sigma(1);

  // compute linearizations w.r.t. lambda_f_sq
  K.Invert();

  const Number dRcauchy_stress_d_lambda_f_sq = evaluate_d_residuum_cauchy_stress_d_lambda_f_sq();
  const Number dRgrowth_scalar_d_lambda_f_sq = evaluate_d_residuum_growth_scalar_d_lambda_f_sq();

  computed_dgrowth_scalar_dlambda_f_sq_ =
      -K(0, 0) * dRgrowth_scalar_d_lambda_f_sq - K(0, 1) * dRcauchy_stress_d_lambda_f_sq;
  computed_dsigma_dlambda_f_sq_ =
      -K(1, 0) * dRgrowth_scalar_d_lambda_f_sq - K(1, 1) * dRcauchy_stress_d_lambda_f_sq;
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::ReinitializeHistory(
    const Number lambda_f, const double time)
{
  ReinitializeState(*this, lambda_f, time);

  const Number dsig = std::invoke(
      [&]()
      {
        if (std::abs(compute_history_cauchy_stress(lambda_f) - sig_h_) < 1e-12)
        {
          return Number(0.0);
        }
        else
        {
          return Number((compute_history_cauchy_stress(lambda_f) - sig_h_) / sig_h_);
        }
      });

  const Number growth_scalar = std::invoke(
      [&]()
      {
        if (history_.size() == 0) return Number(1.0);
        return history_.back().timesteps.back().growth_scalar;
      });

  // initialize history
  const Number growth_scalar_production_rate =
      growth_evolution_.evaluate_true_mass_production_rate(dsig);
  const MassIncrement<Number> mass_increment = {.reference_stretch = lambda_pre_ / lambda_f,
      .growth_scalar = growth_scalar,
      .growth_scalar_production_rate = growth_scalar_production_rate,
      .deposition_time = time};

  if (history_.size() == 0)
  {
    history_.emplace_back();
    history_.back().timesteps.emplace_back(std::move(mass_increment));
    return;
  }

  const MassIncrement<Number> last_mass_increment = history_.back().timesteps.back();
  // check if the item is already in the history
  if (IsAlmostEqual(mass_increment, last_mass_increment, 1e-7))
  {
    return;
  }

  if (!IsNear(mass_increment.deposition_time, last_mass_increment.deposition_time))
  {
    FOUR_C_THROW(
        "The history is not empty, but you reinitialized the fiber with a different deposition "
        "time than the previous one. I don't know what happened in between. This is fatal. You "
        "maybe forgot to adapt the depositon time of all previous times?");
  }

  if (history_.back().timesteps.size() == 1)
  {
    history_.back().timesteps.back() = mass_increment;
    return;
  }

  history_.emplace_back();
  history_.back().timesteps.emplace_back(std::move(mass_increment));

#ifdef FOUR_C_ENABLE_ASSERTIONS
  state_is_set_ = false;
#endif
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::AddTime(const double delta_time)
{
  for (auto& interval : history_)
  {
    for (auto& increment : interval.timesteps)
    {
      increment.deposition_time += delta_time;
    }
  }
  reference_time_ += delta_time;
}

template <typename Number>
double MIXTURE::FullConstrainedMixtureFiber<Number>::get_last_time_in_history() const
{
  if (history_.size() == 0) return reference_time_;
  return history_.back().timesteps.back().deposition_time;
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::Update()
{
  if (growth_enabled_)
  {
    FOUR_C_ASSERT(
        current_time_shift_ == 0.0, "The time shift should be zero if growth is enabled!");
    history_.back().timesteps.emplace_back(
        evaluate_current_mass_increment(computed_growth_scalar_, computed_sigma_));

    // erase depending on some condition
    if (adaptive_history_strategy_ != HistoryAdaptionStrategy::none)
    {
      Number tolerance_per_time =
          adaptive_tolerance_ / (current_time_ - history_[0].timesteps[0].deposition_time);
      switch (adaptive_history_strategy_)
      {
        case HistoryAdaptionStrategy::model_equation:
        {
          for (auto& interval : history_)
          {
            if (interval.timesteps.size() >= 4)
            {
              std::vector<bool> erase_item;
              TimestepAdaptivityInfo adaptivity_info;
              std::tie(erase_item, adaptivity_info) =
                  OptimizeHistoryIntegration(interval.adaptivity_info, interval.timesteps.size(),
                      [&](const std::array<std::optional<unsigned int>, 5>& indices)
                      {
                        const double begin_time = interval.adaptivity_info.GetIndexTime(
                            indices[0].value(), 0.0, interval.base_dt);
                        const double end_time = interval.adaptivity_info.GetIndexTime(
                            indices[4].value(), 0.0, interval.base_dt);
                        return IsModelEquationSimpsonRuleIntegrationBelowTolerance<Number>(
                            growth_evolution_, current_time_, begin_time, end_time,
                            tolerance_per_time * (end_time - begin_time));
                      });

              interval.adaptivity_info = adaptivity_info;


              interval.timesteps.erase(
                  std::remove_if(interval.timesteps.begin(), interval.timesteps.end(),
                      [&](const MassIncrement<Number>& item)
                      { return erase_item.at(&item - interval.timesteps.data()); }),
                  interval.timesteps.end());
            }
          }
          break;
        }
        case HistoryAdaptionStrategy::higher_order_integration:
        {
          for (auto& interval : history_)
          {
            if (interval.timesteps.size() >= 4)
            {
              std::vector<bool> erase_item;
              TimestepAdaptivityInfo adaptivity_info;
              std::tie(erase_item, adaptivity_info) =
                  OptimizeHistoryIntegration(interval.adaptivity_info, interval.timesteps.size(),
                      [&](const std::array<std::optional<unsigned int>, 5>& indices)
                      {
                        // here I need to do a 5th order integration and compare it with 3rd oder
                        auto ComputeIntegrationError = [&](auto integrand)
                        {
                          std::array<std::tuple<double, Number>, 5> values{};
                          for (std::size_t i = 0; i < 5; ++i)
                          {
                            values[i] = std::make_tuple(
                                interval.timesteps[indices[i].value()].deposition_time,
                                integrand(interval.timesteps[indices[i].value()]));
                          }

                          return std::abs(
                              Core::UTILS::IntegrateSimpsonStep(values[0], values[2], values[4]) -
                              IntegrateBooleStep(values));
                        };

                        auto current_growth_scalar_integrand =
                            [&](const MassIncrement<Number>& mass_increment)
                        { return growth_scalar_integrand(mass_increment, current_time_); };

                        auto current_cauchy_stress_integrand =
                            [&](const MassIncrement<Number>& mass_increment) {
                              return scaled_cauchy_stress_integrand(
                                  mass_increment, current_time_, current_state_.lambda_f);
                            };

                        const double begin_time = interval.adaptivity_info.GetIndexTime(
                            indices[0].value(), 0.0, interval.base_dt);
                        const double end_time = interval.adaptivity_info.GetIndexTime(
                            indices[4].value(), 0.0, interval.base_dt);
                        Number allowed_tolerance = tolerance_per_time * (end_time - begin_time);

                        bool coarsen =
                            ComputeIntegrationError(current_growth_scalar_integrand) <=
                                allowed_tolerance &&
                            ComputeIntegrationError(current_cauchy_stress_integrand) / sig_h_ <=
                                allowed_tolerance;
                        return coarsen;
                      });

              interval.adaptivity_info = adaptivity_info;

              interval.timesteps.erase(
                  std::remove_if(interval.timesteps.begin(), interval.timesteps.end(),
                      [&](const MassIncrement<Number>& item)
                      { return erase_item.at(&item - interval.timesteps.data()); }),
                  interval.timesteps.end());
            }
          }
          break;
        }
        case HistoryAdaptionStrategy::window:
        {
          std::size_t num_total_items = 0;
          for (auto i = history_.rbegin(); i != history_.rend(); ++i)
          {
            auto& interval = *i;
            std::vector<bool> erase_item(interval.timesteps.size(), false);

            num_total_items += interval.timesteps.size();
            if (num_total_items > window_size)
            {
              std::size_t num_to_delete = num_total_items - window_size;
              for (std::size_t index = 0;
                   index < std::min(num_to_delete, interval.timesteps.size()); ++index)
              {
                erase_item[index] = true;
              }
            }

            interval.timesteps.erase(
                std::remove_if(interval.timesteps.begin(), interval.timesteps.end(),
                    [&](const MassIncrement<Number>& item)
                    { return erase_item.at(&item - interval.timesteps.data()); }),
                interval.timesteps.end());
          }

          history_.erase(std::remove_if(history_.begin(), history_.end(),
                             [](const auto& item) { return item.timesteps.size() <= 1; }),
              history_.end());
          break;
        }
        default:
          FOUR_C_THROW("unknown history adaption strategy");
      }
    }
  }
  else
  {
    const double delta_time = current_time_ - get_last_time_in_history();
    AddTime(delta_time);
  }

  current_time_shift_ = 0.0;
#ifdef FOUR_C_ENABLE_ASSERTIONS
  state_is_set_ = false;
#endif
}


template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::set_deposition_stretch(double lambda_pre)
{
  lambda_pre_ = lambda_pre;
  sig_h_ = EvaluateFiberMaterialCauchyStress<Number>(*fiber_material_, lambda_pre_);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<
    Number>::evaluate_d_current_fiber_p_k2_stress_d_lambdafsq() const
{
  FOUR_C_ASSERT(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  const Number lambda_f_sq = std::pow(current_state_.lambda_f, 2);
  return (computed_dsigma_dlambda_f_sq_ * lambda_f_sq - computed_sigma_) / std::pow(lambda_f_sq, 2);
}



template class MIXTURE::FullConstrainedMixtureFiber<double>;
template class MIXTURE::FullConstrainedMixtureFiber<Sacado::Fad::DFad<double>>;
FOUR_C_NAMESPACE_CLOSE
