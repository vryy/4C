#ifndef FOUR_C_UTILS_MATHOPERATIONS_CLN_HPP
#define FOUR_C_UTILS_MATHOPERATIONS_CLN_HPP

#include "4C_config.hpp"

#include "4C_utils_clnwrapper.hpp"
#include "4C_utils_mathoperations.hpp"

#include <cln/cln.h>


FOUR_C_NAMESPACE_OPEN

namespace Core
{
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
    static constexpr double get_double(T t) { return cln::double_approx(t.Value()); }
  };

}  // namespace Core
FOUR_C_NAMESPACE_CLOSE

#endif
