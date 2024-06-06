/*----------------------------------------------------------------------*/
/*! \file

\brief Various mathematicl functions that can be applied to quantities of artithmetic type, CLN
type, or Sacado FAD type.

\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_UTILS_MATHOPERATIONS_HPP
#define FOUR_C_UTILS_MATHOPERATIONS_HPP

#include "4C_config.hpp"

#include "4C_utils_clnwrapper.hpp"

#include <cln/cln.h>
#include <Sacado.hpp>

#include <cmath>
#include <cstdlib>

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

  template <typename T>
  struct MathOperations<T, std::enable_if_t<std::is_same_v<std::decay_t<T>, Core::CLN::ClnWrapper>>>
  {
    static constexpr T abs(const T& t) { return cln::abs(t.Value()); }
    static constexpr T sqrt(const T& t) { return cln::sqrt(t.Value()); }
    static constexpr T pow(const T& base, const int exponent)
    {
      FOUR_C_THROW("Not implemented!");
      return base;
    }
    static constexpr double GetDouble(T t) { return cln::double_approx(t.Value()); }
  };

  template <typename T>
  struct MathOperations<T, std::enable_if_t<Sacado::IsFad<std::decay_t<T>>::value>>
  {
    template <typename T2>
    static constexpr std::decay_t<T2> abs(T2&& t)
    {
      // Note: This only works since Sacado exports their math functions to namespace std and not
      // one of their own.
      return std::abs(std::forward<T2>(t));
    }
    template <typename T2>
    static constexpr std::decay_t<T2> sqrt(T2&& t)
    {
      // Note: This only works since Sacado exports their math functions to namespace std and not
      // one of their own.
      return std::sqrt(std::forward<T2>(t));
    }
    template <typename T2>
    static constexpr std::decay_t<T2> pow(T2&& base, int exponent)
    {
      // Note: This only works since Sacado exports their math functions to namespace std and not
      // one of their own.
      return std::pow(std::forward<T2>(base), exponent);
    }
  };
}  // namespace Core
FOUR_C_NAMESPACE_CLOSE

#endif
