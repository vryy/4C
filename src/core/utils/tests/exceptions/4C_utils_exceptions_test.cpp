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

  TEST(ExceptionsTest, Exception)
  {
    int a = 1, b = 0;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        division(a, b), Core::Exception, "Division by zero condition!");
  }

  TEST(ExceptionsTest, AssertAlwaysPrintsTest)
  {
    const auto always = []() { FOUR_C_THROW_UNLESS(1 == 2, "Throw."); };

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(always(), Core::Exception, "1 == 2");
  }
}  // namespace
