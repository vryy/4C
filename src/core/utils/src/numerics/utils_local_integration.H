/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation and helper functions for local integration methods

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef UTILS_LOCAL_INTEGRATION_H
#define UTILS_LOCAL_INTEGRATION_H

#include <functional>
#include <iterator>
#include <tuple>

namespace CORE::UTILS
{
  /*!
   * @brief Integrate the interval of 3 datapoints using the Simpson's rule.
   *
   * @tparam StepData A tuple or pair holding the time and the value of the evaluated function.
   * @param step_data_a Left end of the integration inverval.
   * @param step_data_b Point between left and right end of the integration point
   * @param step_data_c Right point of the integration interval.
   * @return ValueType Result of the integration
   */
  template <typename StepData>
  auto IntegrateSimpsonStep(
      const StepData& step_data_a, const StepData& step_data_b, const StepData& step_data_c)
      -> std::tuple_element_t<1, std::remove_reference_t<decltype(step_data_a)>>
  {
    const auto& [time1, value1] = step_data_a;
    const auto& [time2, value2] = step_data_b;
    const auto& [time3, value3] = step_data_c;

    const auto h0 = time2 - time1;
    const auto h1 = time3 - time2;
    const auto h = time3 - time1;

    const auto alpha = (2 - h1 / h0);
    const auto beta = (h * h / (h0 * h1));
    const auto gamma = (2 - h0 / h1);

    return h / 6.0 * (alpha * value1 + beta * value2 + gamma * value3);
  }

  /*!
   * @brief Integrate the interval between the midpoint and the endpoint using the Simpson's rule,
   * given a third point left to the midpoint.
   *
   * @tparam StepData A type with the tuple<TimeType, ValueType>.
   * @param step_data_a Point left of the beginning of the integration interval.
   * @param step_data_b Left end of the integration inverval.
   * @param step_data_c Right point of the integration interval.
   * @return std::tuple<ValueType, DerivativeType> Result of the integration and the derivative of
   * the integration rule w.r.t. the integrand.
   */
  template <typename StepData>
  static inline auto IntegrateSimpsonStepBCAndReturnDerivativeC(
      const StepData& step_data_a, const StepData& step_data_b, const StepData& step_data_c)
      -> std::tuple<std::tuple_element_t<1, std::remove_reference_t<decltype(step_data_a)>>,
          std::tuple_element_t<0, std::remove_reference_t<decltype(step_data_a)>>>
  {
    using ValueType = std::tuple_element_t<1, std::remove_reference_t<decltype(step_data_a)>>;
    using TimeType = std::tuple_element_t<0, std::remove_reference_t<decltype(step_data_a)>>;
    const auto& [time_a, value_a] = step_data_a;
    const auto& [time_b, value_b] = step_data_b;
    const auto& [time_c, value_c] = step_data_c;

    const auto h0 = time_b - time_a;
    const auto h1 = time_c - time_b;
    const auto h = time_c - time_a;

    const auto h1sq = h1 * h1;
    const auto alpha = (2 * h1sq + 3 * h0 * h1) / (6 * h);
    const auto beta = (h1sq + 3 * h0 * h1) / (6 * h0);
    const auto gamma = -h1sq * h1 / (6 * h0 * h);

    const ValueType value = alpha * value_c + beta * value_b + gamma * value_a;
    const TimeType derivative = alpha;

    return {value, derivative};
  }

  /*!
   * @brief See @p CORE::UTILS::IntegrateSimpsonStepBCAndReturnDerivativeC, without returning the
   * partial derivative.
   *
   * @tparam StepData
   * @param step_data1
   * @param step_data2
   * @param step_data3
   * @return auto
   */
  template <typename StepData>
  auto IntegrateSimpsonStepBC(
      const StepData& step_data_a, const StepData& step_data_b, const StepData& step_data_c)
      -> std::tuple_element_t<1, std::remove_reference_t<decltype(step_data_a)>>
  {
    return std::get<0>(
        IntegrateSimpsonStepBCAndReturnDerivativeC(step_data_a, step_data_b, step_data_c));
  }

  /*!
   * @brief Integrate the interval between the startpoint and the endpoint using the Trapezoidal
   * rule.
   *
   * @tparam StepData A tuple containing the time and the value of the integrand
   * @param step_data_a Left end of the integration inverval.
   * @param step_data_b Right end of the integration inverval.
   * @return A function object that takes @p StepData and returns a tuple of the time and
   * the value.
   * @return std::tuple<ValueType, DerivativeType> Result of the integration and the derivative of
   * the integration rule w.r.t. the integrand.
   */
  template <typename StepData>
  auto IntegrateTrapezoidalStepAndReturnDerivativeB(
      const StepData& step_data_a, const StepData& step_data_b)
      -> std::tuple<std::tuple_element_t<1, std::remove_reference_t<decltype(step_data_a)>>,
          std::tuple_element_t<0, std::remove_reference_t<decltype(step_data_a)>>>
  {
    using ValueType = std::tuple_element_t<1, std::remove_reference_t<decltype(step_data_a)>>;
    using TimeType = std::tuple_element_t<0, std::remove_reference_t<decltype(step_data_a)>>;
    const auto& [time_a, value_a] = step_data_a;
    const auto& [time_b, value_b] = step_data_b;


    const ValueType value = 0.5 * (time_b - time_a) * (value_a + value_b);
    const TimeType derivative = 0.5 * (time_b - time_a);

    return std::make_tuple(value, derivative);
  }

  /*!
   * @brief See @p CORE::UTILS::IntegrateTrapezoidalStepAndReturnDerivativeB, without returning the
   * partial derivative.
   *
   * @tparam StepData
   * @param step_data_a
   * @param step_data_b
   * @return auto
   */
  template <typename StepData>
  auto IntegrateTrapezoidalStep(const StepData& step_data_a, const StepData& step_data_b)
      -> std::tuple_element_t<1, std::remove_reference_t<decltype(step_data_a)>>
  {
    return std::get<0>(IntegrateTrapezoidalStepAndReturnDerivativeB(step_data_a, step_data_b));
  }

  /*!
   * @brief Integrate over the given snapshots and the integrator using the Trapezoidal Rule if the
   * size is 2 or the Simpson's rule if the size is larger or equal 3.
   *
   * @tparam Integrand
   * @tparam StepDataIt
   * @param indexable_container An indexable container (e.g. std::vector)
   * @param integrand The integrand evaluated at the elements
   * @return ValueType Integration result
   */
  template <typename Integrand, typename Container>
  inline auto IntegrateSimpsonTrapezoidal(Container indexable_container, Integrand integrand)
      -> std::tuple_element_t<1, decltype(integrand(indexable_container[0]))>
  {
    using TupleType = decltype(integrand(indexable_container[0]));
    using ValueType = std::tuple_element_t<1, TupleType>;
    auto size = std::size(indexable_container);

    if (size <= 1)
    {
      // nothing to integrate
      return 0;
    }

    if (size == 2)
    {
      // Apply trapezoidal rule
      return IntegrateTrapezoidalStep<TupleType>(
          integrand(indexable_container[0]), integrand(indexable_container[1]));
    }

    // apply Simpson's rule
    ValueType integration_result = 0;

    TupleType last_evaluation = integrand(indexable_container[0]);
    // Apply simpson's rule
    for (unsigned i = 1; i < size - 1; ++ ++i)
    {
      const TupleType mid_evaluation = integrand(indexable_container[i]);
      const TupleType end_evaluation = integrand(indexable_container[i + 1]);
      integration_result +=
          IntegrateSimpsonStep<TupleType>(last_evaluation, mid_evaluation, end_evaluation);
      last_evaluation = end_evaluation;
    }

    // special case for even number
    if (size % 2 == 0)
    {
      integration_result +=
          IntegrateSimpsonStepBC<TupleType>(integrand(indexable_container[size - 3]),
              last_evaluation, integrand(indexable_container[size - 1]));
    }

    return integration_result;
  }
}  // namespace CORE::UTILS

#endif
