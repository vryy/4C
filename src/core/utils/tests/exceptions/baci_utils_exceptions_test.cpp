/*----------------------------------------------------------------------*/
/*! \file
 *
\brief Testcases for the core utils

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_utils_exceptions.hpp"

#include "baci_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace FourC;

  int division(int a, int b)
  {
    if (b == 0)
    {
      FOUR_C_THROW("Division by zero condition!");
    }
    return (a / b);
  }

  TEST(CoreUtilsTest, Exception)
  {
    int a = 1, b = 0;

    BACI_EXPECT_THROW_WITH_MESSAGE(division(a, b), CORE::Exception, "Division by zero condition!");
  }
}  // namespace
