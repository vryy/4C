/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a 1D full constrained mixture fiber
\level 3
*/
/*----------------------------------------------------------------------*/

#include "mixture_full_constrained_mixture_fiber.H"
#include <algorithm>
#include <limits>
#include <numeric>
#include <memory>
#include "utils_exceptions.H"
#include "lib_pack_buffer.H"
#include "mixture_constituent_remodelfiber_material.H"
#include "lib_parobject.H"
#include <Sacado.hpp>
#include <type_traits>
#include "utils_local_integration.H"
#include "utils_local_newton.H"
#include "mixture_growth_evolution_linear_cauchy_poisson_turnover.H"

namespace
{
  template <typename Number>
  static inline Number ComputeAbsoluteErrorWhenSkippingSnapshot(
      std::size_t i, const std::vector<std::tuple<double, Number>>& evaluated_integrand)
  {
    dsassert(i > 0,
        "This function should only be called for an item in the interior such that it can still be "
        "integrated with a Simpson's rule if we remove point i (0 < i < size-2)");
    dsassert(i < evaluated_integrand.size() - 2,
        "This function should only be called for an item in the interior such that it can still be "
        "integrated with a Simpson's rule if we remove point i (0 < i < size-2)");

    Number full_integration =
        CORE::UTILS::IntegrateSimpsonStep<std::tuple<double, Number>>(
            evaluated_integrand[i - 1], evaluated_integrand[i], evaluated_integrand[i + 1]) +
        CORE::UTILS::IntegrateSimpsonStepBC<std::tuple<double, Number>>(
            evaluated_integrand[i + 0], evaluated_integrand[i + 1], evaluated_integrand[i + 2]);
    Number skipped_integration = CORE::UTILS::IntegrateSimpsonStep<std::tuple<double, Number>>(
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
      integration_result += CORE::UTILS::IntegrateSimpsonTrapezoidal(interval,
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
    const size_t size = interval.size();
    if (size == 0) dserror("The history is empty. I cannot integrate.");
    if (size == 1)
    {
      // can only apply trapezoidal rule
      const auto [integration, derivative] =
          CORE::UTILS::IntegrateTrapezoidalStepAndReturnDerivativeB<std::tuple<double, Number>>(
              {interval[0].deposition_time, integrand(interval[0])},
              {current_increment.deposition_time, integrand(current_increment)});
      return std::make_tuple(integration + history_integration, derivative);
    }

    const auto [integration, derivative] =
        CORE::UTILS::IntegrateSimpsonStepBCAndReturnDerivativeC<std::tuple<double, Number>>(
            {interval[size - 2].deposition_time, integrand(interval[size - 2])},
            {interval[size - 1].deposition_time, integrand(interval[size - 1])},
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
#ifdef DEBUG
    fiber.state_is_set_ = true;
#endif
    fiber.current_state_.lambda_f = lambda_f;
    fiber.current_time_ = time;
  }
}  // namespace

template <typename Number>
MIXTURE::FullConstrainedMixtureFiber<Number>::FullConstrainedMixtureFiber(
    std::shared_ptr<const MIXTURE::RemodelFiberMaterial<Number>> material,
    MIXTURE::LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<Number> growth_evolution,
    Number lambda_pre, bool adaptive_history, bool growth_enabled)
    : lambda_pre_(lambda_pre),
      fiber_material_(std::move(material)),
      growth_evolution_(growth_evolution),
      growth_enabled_(growth_enabled),
      adaptive_history_(adaptive_history),
      current_time_(0),
      computed_growth_scalar_(1.0)
{
  sig_h_ = EvaluateFiberMaterialCauchyStress<Number>(*fiber_material_, lambda_pre_);
  computed_sigma_ = sig_h_;
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::Pack(DRT::PackBuffer& data) const
{
  dserror("Packing and Unpacking is currently only implemented for the double-specialization");
}

template <>
void MIXTURE::FullConstrainedMixtureFiber<double>::Pack(DRT::PackBuffer& data) const
{
  data.AddtoPack(sig_h_);
  data.AddtoPack(lambda_pre_);

  data.AddtoPack(current_state_.lambda_f);

  data.AddtoPack(reference_time_);
  data.AddtoPack(current_time_shift_);

  data.AddtoPack(history_.size());
  for (const auto& interval : history_)
  {
    data.AddtoPack(interval.size());
    for (const auto& item : interval)
    {
      data.AddtoPack(item.reference_stretch);
      data.AddtoPack(item.growth_scalar);
      data.AddtoPack(item.growth_scalar_production_rate);
      data.AddtoPack(item.deposition_time);
    }
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
  dserror("Packing and Unpacking is currently only implemented for the double-specialization");
}

template <>
void MIXTURE::FullConstrainedMixtureFiber<double>::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  DRT::ParObject::ExtractfromPack(position, data, sig_h_);
  DRT::ParObject::ExtractfromPack(position, data, lambda_pre_);
  DRT::ParObject::ExtractfromPack(position, data, current_state_.lambda_f);

  DRT::ParObject::ExtractfromPack(position, data, reference_time_);
  DRT::ParObject::ExtractfromPack(position, data, current_time_shift_);

  std::size_t size_of_history;
  DRT::ParObject::ExtractfromPack(position, data, size_of_history);
  history_.resize(size_of_history);

  for (auto& interval : history_)
  {
    std::size_t size_of_interval;
    DRT::ParObject::ExtractfromPack(position, data, size_of_interval);
    interval.resize(size_of_interval);
    for (auto& item : interval)
    {
      DRT::ParObject::ExtractfromPack(position, data, item.reference_stretch);
      DRT::ParObject::ExtractfromPack(position, data, item.growth_scalar);
      DRT::ParObject::ExtractfromPack(position, data, item.growth_scalar_production_rate);
      DRT::ParObject::ExtractfromPack(position, data, item.deposition_time);
    }
  }


  DRT::ParObject::ExtractfromPack(position, data, current_time_);


  DRT::ParObject::ExtractfromPack(position, data, computed_growth_scalar_);
  DRT::ParObject::ExtractfromPack(position, data, computed_sigma_);
  DRT::ParObject::ExtractfromPack(position, data, computed_dgrowth_scalar_dlambda_f_sq_);
  DRT::ParObject::ExtractfromPack(position, data, computed_dsigma_dlambda_f_sq_);
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::RecomputeState(
    const Number lambda_f, const double time)
{
  ReinitializeState(*this, lambda_f, time);

  if (!growth_enabled_)
  {
    current_time_shift_ = GetLastTimeInHistory() - current_time_;
  }
  else
  {
    current_time_shift_ = 0.0;
  }

  ComputeInternalVariables();
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::ComputeHistoryCauchyStress(
    const Number lambda_f) const
{
  const Number growth_scalar = std::invoke(
      [&]()
      {
        if (history_.size() > 0) return history_.back().back().growth_scalar;
        return Number(1.0);
      });
  const double time = std::invoke(
      [&]()
      {
        if (history_.size() > 0) return history_.back().back().deposition_time;
        return 0.0;
      });


  const auto scaled_cauchy_stress_integrand = [&](const MassIncrement<Number>& mass_increment)
  { return ScaledCauchyStressIntegrand(mass_increment, time, lambda_f); };

  const Number scaled_cauchy_stress_history =
      growth_evolution_.EvaluateSurvivalFunction(time) *
          EvaluateFiberMaterialCauchyStress<Number>(*fiber_material_, lambda_pre_ * lambda_f) +
      IntegrateOverDepositionHistory<Number>(history_, scaled_cauchy_stress_integrand);

  return scaled_cauchy_stress_history / growth_scalar;
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::GrowthScalarIntegrand(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time) const
{
  return mass_increment.growth_scalar_production_rate * mass_increment.growth_scalar *
         growth_evolution_.EvaluateSurvivalFunction(time - mass_increment.deposition_time);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::DGrowthScalarIntegrandDProductionRate(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time) const
{
  return mass_increment.growth_scalar *
         growth_evolution_.EvaluateSurvivalFunction(time - mass_increment.deposition_time);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::DGrowthScalarIntegrandDGrowthScalar(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time) const
{
  return mass_increment.growth_scalar_production_rate *
         growth_evolution_.EvaluateSurvivalFunction(time - mass_increment.deposition_time);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::ScaledCauchyStressIntegrand(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  return mass_increment.growth_scalar_production_rate * mass_increment.growth_scalar *
         growth_evolution_.EvaluateSurvivalFunction(time - mass_increment.deposition_time) *
         EvaluateFiberMaterialCauchyStress<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::DScaledCauchyStressIntegrandDProductionRate(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  dsassert(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  return mass_increment.growth_scalar *
         growth_evolution_.EvaluateSurvivalFunction(time - mass_increment.deposition_time) *
         EvaluateFiberMaterialCauchyStress<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::DScaledCauchyStressIntegrandDGrowthScalar(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  dsassert(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  return mass_increment.growth_scalar_production_rate *
         growth_evolution_.EvaluateSurvivalFunction(time - mass_increment.deposition_time) *
         EvaluateFiberMaterialCauchyStress<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::DScaledCauchyStressIntegrandDLambdaFSq(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  dsassert(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  return mass_increment.growth_scalar_production_rate * mass_increment.growth_scalar *
         growth_evolution_.EvaluateSurvivalFunction(time - mass_increment.deposition_time) *
         EvaluateDFiberMaterialCauchyStressDLambdaESq<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f) *
         mass_increment.reference_stretch * mass_increment.reference_stretch;
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::DScaledCauchyStressIntegrandDLambdaRefSq(
    const MIXTURE::MassIncrement<Number>& mass_increment, const double time,
    const Number current_lambda_f) const
{
  dsassert(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  return mass_increment.growth_scalar_production_rate * mass_increment.growth_scalar *
         growth_evolution_.EvaluateSurvivalFunction(time - mass_increment.deposition_time) *
         EvaluateDFiberMaterialCauchyStressDLambdaESq<Number>(
             *fiber_material_, mass_increment.reference_stretch * current_lambda_f) *
         current_lambda_f * current_lambda_f;
}

template <typename Number>
std::function<std::tuple<LINALG::Matrix<2, 1, Number>, LINALG::Matrix<2, 2, Number>>(
    const LINALG::Matrix<2, 1, Number>&)>
MIXTURE::FullConstrainedMixtureFiber<Number>::GetLocalNewtonEvaluator() const
{
  dsassert(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  dsassert(current_time_shift_ == 0.0, "The timeshift should be zero if growth is enabled");
  dsassert(
      history_.size() > 0, "The history is empty. Please initialize it with ReinitializeHistory()");

  const auto growth_scalar_integrand = [&](const MassIncrement<Number>& mass_increment)
  { return GrowthScalarIntegrand(mass_increment, current_time_); };

  const auto dgrowth_scalar_integrand_dproduction_rate =
      [&](const MassIncrement<Number>& mass_increment)
  { return DGrowthScalarIntegrandDProductionRate(mass_increment, current_time_); };

  const auto dgrowth_scalar_integrand_dgrowth_scalar =
      [&](const MassIncrement<Number>& mass_increment)
  { return DGrowthScalarIntegrandDGrowthScalar(mass_increment, current_time_); };

  const auto scaled_cauchy_stress_integrand = [&](const MassIncrement<Number>& mass_increment)
  { return ScaledCauchyStressIntegrand(mass_increment, current_time_, current_state_.lambda_f); };

  const auto dscaled_cauchy_stress_integrand_dproduction_rate =
      [&](const MassIncrement<Number>& mass_increment)
  {
    return DScaledCauchyStressIntegrandDProductionRate(
        mass_increment, current_time_, current_state_.lambda_f);
  };

  const auto dscaled_cauchy_stress_integrand_dgrowth_scalar =
      [&](const MassIncrement<Number>& mass_increment)
  {
    return DScaledCauchyStressIntegrandDGrowthScalar(
        mass_increment, current_time_, current_state_.lambda_f);
  };

  // evaluate history quantities (except last timestep, which depends on the current Cauchy
  // stress)
  const Number growth_scalar_history =
      growth_evolution_.EvaluateSurvivalFunction(current_time_ - reference_time_) +
      IntegrateOverDepositionHistory<Number>(history_, growth_scalar_integrand);

  const Number scaled_cauchy_stress_history =
      growth_evolution_.EvaluateSurvivalFunction(current_time_ - reference_time_) *
          EvaluateFiberMaterialCauchyStress<Number>(
              *fiber_material_, lambda_pre_ * current_state_.lambda_f) +
      IntegrateOverDepositionHistory<Number>(history_, scaled_cauchy_stress_integrand);

  return [=](const LINALG::Matrix<2, 1, Number>& growth_scalar_and_cauchy_stress)
  {
    const Number growth_scalar = growth_scalar_and_cauchy_stress(0);
    const Number cauchy_stress = growth_scalar_and_cauchy_stress(1);

    const MassIncrement<Number> current_increment =
        EvaluateCurrentMassIncrement(growth_scalar, cauchy_stress);

    const Number dmy_growth_scalar_production_rate_dsig =
        growth_evolution_.EvaluateDTrueMassProductionRateDSig((cauchy_stress - sig_h_) / sig_h_) /
        sig_h_;
    const auto [my_growth_scalar, dmy_growth_scalar_dintegrand] =
        IntegrateLastTimestepWithDerivative<Number>(
            history_.back(), current_increment, growth_scalar_integrand, growth_scalar_history);

    const Number dmy_growth_scalar_dsig =
        dmy_growth_scalar_dintegrand *
        dgrowth_scalar_integrand_dproduction_rate(current_increment) *
        dmy_growth_scalar_production_rate_dsig;
    const Number dmy_growth_scalar_dgrowth_scalar =
        dmy_growth_scalar_dintegrand * dgrowth_scalar_integrand_dgrowth_scalar(current_increment);

    const auto [my_scaled_cauchy_stress, dmy_scaled_cauchy_stress_dintegrand] =
        IntegrateLastTimestepWithDerivative<Number>(history_.back(), current_increment,
            scaled_cauchy_stress_integrand, scaled_cauchy_stress_history);


    const Number dmy_scaled_cauchy_stress_dsig =
        dmy_scaled_cauchy_stress_dintegrand *
        dscaled_cauchy_stress_integrand_dproduction_rate(current_increment) *
        dmy_growth_scalar_production_rate_dsig;
    const Number dmy_scaled_cauchy_stress_dgrowth_scalar =
        dmy_scaled_cauchy_stress_dintegrand *
        dscaled_cauchy_stress_integrand_dgrowth_scalar(current_increment);

    LINALG::Matrix<2, 1, Number> residua;
    residua(0) = my_growth_scalar - growth_scalar;
    residua(1) = my_scaled_cauchy_stress / growth_scalar - cauchy_stress;

    LINALG::Matrix<2, 2, Number> derivative;
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
Number MIXTURE::FullConstrainedMixtureFiber<Number>::EvaluateDResiduumGrowthScalarDLambdaFSq() const
{
  return 0.0;
}


template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::EvaluateDResiduumCauchyStressDLambdaFSq() const
{
  dsassert(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  dsassert(current_time_shift_ == 0.0, "The timeshift should be zero if growth is enabled");
  dsassert(
      history_.size() > 0, "The history is empty. Please initialize it with ReinitializeHistory()");

  const auto dscaled_cauchy_stress_integrand_d_lambda_f_sq =
      [&](const MassIncrement<Number>& mass_increment)
  {
    return DScaledCauchyStressIntegrandDLambdaFSq(
        mass_increment, current_time_, current_state_.lambda_f);
  };
  const Number Dscaled_cauchy_stress_D_lambda_f_sq_history =
      growth_evolution_.EvaluateSurvivalFunction(current_time_ - reference_time_) *
          EvaluateDFiberMaterialCauchyStressDLambdaESq<Number>(
              *fiber_material_, lambda_pre_ * current_state_.lambda_f) *
          std::pow(lambda_pre_, 2) +
      IntegrateOverDepositionHistory<Number>(
          history_, dscaled_cauchy_stress_integrand_d_lambda_f_sq);

  const MassIncrement<Number> current_increment =
      EvaluateCurrentMassIncrement(computed_growth_scalar_, computed_sigma_);

  const auto [dscaled_cauchy_stress_lambda_f_sq,
      ddscaled_cauchy_stress_integrand_d_lambda_f_sq_d_integrand] =
      IntegrateLastTimestepWithDerivative(history_.back(), current_increment,
          dscaled_cauchy_stress_integrand_d_lambda_f_sq,
          Dscaled_cauchy_stress_D_lambda_f_sq_history);


  return (dscaled_cauchy_stress_lambda_f_sq +
             ddscaled_cauchy_stress_integrand_d_lambda_f_sq_d_integrand *
                 DScaledCauchyStressIntegrandDLambdaRefSq(
                     current_increment, current_time_, current_state_.lambda_f) *
                 EvaluateDLambdaRefSqDLambdaFSq(current_state_.lambda_f)) /
         computed_growth_scalar_;
}


template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::EvaluateLambdaRef(Number lambda_f) const
{
  return lambda_pre_ / lambda_f;
}


template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::EvaluateDLambdaRefSqDLambdaFSq(
    Number lambda_f) const
{
  return -std::pow(lambda_pre_, 2) / std::pow(lambda_f, 4);
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::ComputeInternalVariables()
{
  if (!growth_enabled_)
  {
    // in this case, I don't need to solve a local Newton method. I just need to do the integration
    // over the history until the moment.
    computed_dgrowth_scalar_dlambda_f_sq_ = 0;

    const auto growth_scalar_integrand = [&](const MassIncrement<Number>& mass_increment)
    { return GrowthScalarIntegrand(mass_increment, current_time_ + current_time_shift_); };

    const auto scaled_cauchy_stress_integrand = [&](const MassIncrement<Number>& mass_increment)
    {
      return ScaledCauchyStressIntegrand(
          mass_increment, current_time_ + current_time_shift_, current_state_.lambda_f);
    };

    const auto dscaled_cauchy_stress_integrand_d_lambda_f_sq =
        [&](const MassIncrement<Number>& mass_increment)
    {
      return DScaledCauchyStressIntegrandDLambdaFSq(
          mass_increment, current_time_, current_state_.lambda_f);
    };

    computed_growth_scalar_ =
        growth_evolution_.EvaluateSurvivalFunction(
            current_time_ + current_time_shift_ - reference_time_) +
        IntegrateOverDepositionHistory<Number>(history_, growth_scalar_integrand);

    computed_sigma_ =
        (growth_evolution_.EvaluateSurvivalFunction(
             current_time_ + current_time_shift_ - reference_time_) *
                EvaluateFiberMaterialCauchyStress<Number>(
                    *fiber_material_, lambda_pre_ * current_state_.lambda_f) +
            IntegrateOverDepositionHistory<Number>(history_, scaled_cauchy_stress_integrand)) /
        computed_growth_scalar_;


    computed_dsigma_dlambda_f_sq_ =
        (growth_evolution_.EvaluateSurvivalFunction(
             current_time_ + current_time_shift_ - reference_time_) *
                EvaluateDFiberMaterialCauchyStressDLambdaESq<Number>(
                    *fiber_material_, lambda_pre_ * current_state_.lambda_f) *
                std::pow(lambda_pre_, 2) +
            IntegrateOverDepositionHistory<Number>(
                history_, dscaled_cauchy_stress_integrand_d_lambda_f_sq)) /
        computed_growth_scalar_;

    return;
  }

  //  Evaluate local Newton system to obtain the current Cauchy stress
  const auto EvaluateCurrentLocalNewtonLinearSystem = GetLocalNewtonEvaluator();

  constexpr auto tolerance = 1e-10;
  constexpr auto max_iterations = 500;

  LINALG::Matrix<2, 1, Number> initial_guess;
  initial_guess(0) = computed_growth_scalar_;
  initial_guess(1) = computed_sigma_;
  auto [growth_scalar_and_sigma, K] = CORE::UTILS::SolveLocalNewtonAndReturnJacobian(
      EvaluateCurrentLocalNewtonLinearSystem, initial_guess, tolerance, max_iterations);

  computed_growth_scalar_ = growth_scalar_and_sigma(0);
  computed_sigma_ = growth_scalar_and_sigma(1);

  // compute linearizations w.r.t. lambda_f_sq
  K.Invert();

  const Number dRcauchy_stress_d_lambda_f_sq = EvaluateDResiduumCauchyStressDLambdaFSq();
  const Number dRgrowth_scalar_d_lambda_f_sq = EvaluateDResiduumGrowthScalarDLambdaFSq();

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
        if (std::abs(ComputeHistoryCauchyStress(lambda_f) - sig_h_) < 1e-12)
        {
          return Number(0.0);
        }
        else
        {
          return Number((ComputeHistoryCauchyStress(lambda_f) - sig_h_) / sig_h_);
        }
      });

  const Number growth_scalar = std::invoke(
      [&]()
      {
        if (history_.size() == 0) return Number(1.0);
        return history_.back().back().growth_scalar;
      });

  // initialize history
  const Number growth_scalar_production_rate =
      growth_evolution_.EvaluateTrueMassProductionRate(dsig);
  const MassIncrement<Number> mass_increment = {.reference_stretch = lambda_pre_ / lambda_f,
      .growth_scalar = growth_scalar,
      .growth_scalar_production_rate = growth_scalar_production_rate,
      .deposition_time = time};

  if (history_.size() == 0)
  {
    history_.emplace_back(std::vector<MassIncrement<Number>>({std::move(mass_increment)}));
    return;
  }

  const MassIncrement<Number> last_mass_increment = history_.back().back();
  // check if the history is already
  if (IsAlmostEqual(mass_increment, last_mass_increment, 1e-7))
  {
    return;
  }

  if (std::abs(mass_increment.deposition_time - last_mass_increment.deposition_time) > 1e-7)
  {
    dserror(
        "The history is not empty, but you reinitialized the fiber with a different deposition "
        "time than the previous one. I don't know what happened in between. This is fatal. You "
        "maybe forgot to adapt the depositon time of all previous times?");
  }

  if (history_.back().size() == 1)
  {
    history_.back().back() = mass_increment;
    return;
  }

  history_.emplace_back(std::vector<MassIncrement<Number>>({mass_increment}));

#ifdef DEBUG
  state_is_set_ = false;
#endif
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::AddTime(const double delta_time)
{
  for (auto& interval : history_)
  {
    for (auto& increment : interval)
    {
      increment.deposition_time += delta_time;
    }
  }
  reference_time_ += delta_time;
}

template <typename Number>
double MIXTURE::FullConstrainedMixtureFiber<Number>::GetLastTimeInHistory() const
{
  if (history_.size() == 0) return reference_time_;
  return history_.back().back().deposition_time;
}

template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::Update()
{
  if (growth_enabled_)
  {
    dsassert(current_time_shift_ == 0.0, "The time shift should be zero if growth is enabled!");
    history_.back().emplace_back(
        EvaluateCurrentMassIncrement(computed_growth_scalar_, computed_sigma_));

    // erase depending on some condition
    if (adaptive_history_)
    {
      for (auto& interval : history_)
      {
        if (interval.size() >= 4)
        {
          std::vector<std::tuple<double, Number>> evaluated_scaled_cauchy_stress_integrand(
              interval.size());
          std::transform(interval.begin(), interval.end(),
              evaluated_scaled_cauchy_stress_integrand.begin(),
              [&](const MassIncrement<Number>& inc)
              {
                return std::make_tuple(inc.deposition_time,
                    ScaledCauchyStressIntegrand(inc, current_time_, current_state_.lambda_f));
              });

          std::vector<std::tuple<double, Number>> evaluatd_growth_scalar_integrand(interval.size());
          std::transform(interval.begin(), interval.end(), evaluatd_growth_scalar_integrand.begin(),
              [&](const MassIncrement<Number>& inc) {
                return std::make_tuple(
                    inc.deposition_time, GrowthScalarIntegrand(inc, current_time_));
              });

          std::vector<Number> relative_error_sum(interval.size(), 0);
          relative_error_sum[0] = std::numeric_limits<double>::max();
          relative_error_sum[interval.size() - 1] = std::numeric_limits<double>::max();
          relative_error_sum[interval.size() - 2] = std::numeric_limits<double>::max();
          for (std::size_t i = 1; i < interval.size() - 2; ++i)
          {
            // compute for every point in the history the error of each integral if we skip this
            // point
            relative_error_sum[i] = ComputeAbsoluteErrorWhenSkippingSnapshot(
                                        i, evaluated_scaled_cauchy_stress_integrand) /
                                    computed_growth_scalar_ / computed_sigma_;
            relative_error_sum[i] +=
                ComputeAbsoluteErrorWhenSkippingSnapshot(i, evaluatd_growth_scalar_integrand) /
                computed_growth_scalar_;
          }

          std::vector<bool> erase_item =
              GetEraseItemsBasedOnRelativeError(relative_error_sum, adaptive_tolerance_);
          interval.erase(std::remove_if(interval.begin(), interval.end(),
                             [&](const MassIncrement<Number>& item)
                             { return erase_item.at(&item - interval.data()); }),
              interval.end());
        }
      }
    }
  }
  else
  {
    const double delta_time = current_time_ - GetLastTimeInHistory();
    AddTime(delta_time);
  }

  current_time_shift_ = 0.0;
#ifdef DEBUG
  state_is_set_ = false;
#endif
}


template <typename Number>
void MIXTURE::FullConstrainedMixtureFiber<Number>::SetDepositionStretch(double lambda_pre)
{
  lambda_pre_ = lambda_pre;
  sig_h_ = EvaluateFiberMaterialCauchyStress<Number>(*fiber_material_, lambda_pre_);
}

template <typename Number>
Number MIXTURE::FullConstrainedMixtureFiber<Number>::EvaluateDCurrentFiberPK2StressDLambdafsq()
    const
{
  dsassert(state_is_set_, "You need to call RecomputeState(...) before calling this method!");
  const Number lambda_f_sq = std::pow(current_state_.lambda_f, 2);
  return (computed_dsigma_dlambda_f_sq_ * lambda_f_sq - computed_sigma_) / std::pow(lambda_f_sq, 2);
}



template class MIXTURE::FullConstrainedMixtureFiber<double>;
template class MIXTURE::FullConstrainedMixtureFiber<Sacado::Fad::DFad<double>>;