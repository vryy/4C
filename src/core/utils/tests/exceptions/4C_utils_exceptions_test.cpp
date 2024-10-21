// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_utils_exceptions.hpp"

#include "4C_unittest_utils_assertions_test.hpp"

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

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        division(a, b), Core::Exception, "Division by zero condition!");
  }
}  // namespace
