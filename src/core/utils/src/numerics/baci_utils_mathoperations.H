/*----------------------------------------------------------------------*/
/*! \file

\brief Various mathematicl functions that can be applied to quantities of artithmetic type, CLN
type, or Sacado FAD type.

\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_UTILS_MATHOPERATIONS_H
#define BACI_UTILS_MATHOPERATIONS_H

#include "baci_config.H"

#include "baci_utils_clnwrapper.H"

#include <cln/cln.h>
#include <Sacado.hpp>

#include <cmath>
#include <cstdlib>


BACI_NAMESPACE_OPEN
namespace CORE
{
  template <typename NumberType, typename Enable = void>
  struct MathOperations
  {
    static constexpr NumberType abs(NumberType t) = delete;
  };

  template <typename T>
  struct MathOperations<T, std::enable_if_t<std::is_arithmetic_v<T>>>
  {
    static constexpr T abs(T t) { return std::abs(t); }
    static constexpr T sqrt(T t) { return std::sqrt(t); }
    static constexpr T pow(T base, int exponent) { return std::pow(base, exponent); }
  };

  template <typename T>
  struct MathOperations<T, std::enable_if_t<std::is_same_v<std::decay_t<T>, CORE::CLN::ClnWrapper>>>
  {
    static constexpr T abs(const T& t) { return cln::abs(t.Value()); }
    static constexpr T sqrt(const T& t) { return cln::sqrt(t.Value()); }
    static constexpr T pow(const T& base, const int exponent) { return cln::pow(base, exponent); }
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
}  // namespace CORE
BACI_NAMESPACE_CLOSE

#endif
