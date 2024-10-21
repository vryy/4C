#ifndef FOUR_C_UTILS_MATHOPERATIONS_HPP
#define FOUR_C_UTILS_MATHOPERATIONS_HPP

#include "4C_config.hpp"

#include <cmath>
#include <cstdlib>
#include <type_traits>

FOUR_C_NAMESPACE_OPEN
namespace Core
{
  /**
   * This struct allows to evaluate basic mathematical functions on values of @p NumberType. The
   * implementations are within the specializations for the different @p NumberType.
   * @tparam Enable This template parameter exists to allow SFINAE
   */
  template <typename NumberType, typename Enable = void>
  struct MathOperations
  {
    // Defer evaluation of static_assert until user tries an instantiation.
    static_assert(std::is_same_v<NumberType, void>,
        "No MathOperations specialization for this type! You need to include the header with the "
        "MathOperations in the translation unit where you encounter this error.");
    //! The absolute value of @p t
    static constexpr NumberType abs(NumberType t) = delete;
    //! The square root value of @p t
    static constexpr NumberType sqrt(NumberType t) = delete;
    //! The power of @p
    static constexpr NumberType pow(const NumberType& base, const int exponent) = delete;
  };

  template <typename T>
  struct MathOperations<T, std::enable_if_t<std::is_arithmetic_v<T>>>
  {
    static constexpr T abs(T t) { return std::abs(t); }
    static constexpr T sqrt(T t) { return std::sqrt(t); }
    static constexpr T pow(T base, int exponent) { return std::pow(base, exponent); }
  };

}  // namespace Core
FOUR_C_NAMESPACE_CLOSE

#endif
