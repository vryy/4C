// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRIC_SEARCH_UTILS_TEST_HPP
#define FOUR_C_GEOMETRIC_SEARCH_UTILS_TEST_HPP

#include <gtest/gtest.h>

#include <cfenv>

namespace
{
  class TestWithFloatingPointExceptions : public ::testing::Test
  {
   public:
    TestWithFloatingPointExceptions() = default;

   protected:
    void SetUp() override
    {
      // Start every test with a clean floating-point environment.
      std::feclearexcept(FE_ALL_EXCEPT);
    }

    void TearDown() override
    {
      // We don't allow a call to ArborX to raise any of the following floating point exceptions.
      constexpr int critical_exceptions = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW;
      const int exceptions = std::fetestexcept(critical_exceptions);
      EXPECT_EQ(exceptions & FE_INVALID, 0) << "FE_INVALID was raised";
      EXPECT_EQ(exceptions & FE_DIVBYZERO, 0) << "FE_DIVBYZERO was raised";
      EXPECT_EQ(exceptions & FE_OVERFLOW, 0) << "FE_OVERFLOW was raised";
    }
  };
}  // namespace

#endif
