// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_utils_std_cxx20_ranges.hpp"

#include <list>

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <typename Container>
  class StdCxx20RangesViewsAll : public testing::Test
  {
   protected:
    [[nodiscard]] Container fill() const { return Container{1, 2, 3, 4, 5}; }
  };


  using TestTypes = ::testing::Types<std::array<int, 5>, std::vector<int>, std::list<int>>;
  TYPED_TEST_SUITE(StdCxx20RangesViewsAll, TestTypes);

  TYPED_TEST(StdCxx20RangesViewsAll, testIterator)
  {
    TypeParam container = this->fill();
    int counter = 0;

    for (const auto element : std_20::ranges::views::all(container))
    {
      (void)element;
      counter++;
    }

    EXPECT_EQ(counter, 5);
  }

  TYPED_TEST(StdCxx20RangesViewsAll, testEleAccess)
  {
    const TypeParam container = this->fill();
    int sum_of_elements = 0;

    for (const auto element : std_20::ranges::views::all(container))
    {
      sum_of_elements += element;
    }

    EXPECT_EQ(sum_of_elements, 15);
  }

  TYPED_TEST(StdCxx20RangesViewsAll, testEleManipulation)
  {
    TypeParam container = this->fill();

    for (auto &element : std_20::ranges::views::all(container))
    {
      element += element;
    }

    EXPECT_EQ(container, (TypeParam{2, 4, 6, 8, 10}));
  }


  template <typename Container>
  class StdCxx20RangesViewsFilter : public testing::Test
  {
   protected:
    [[nodiscard]] Container fill() const { return Container{1, 2, 3, 4, 5}; }
  };


  TYPED_TEST_SUITE(StdCxx20RangesViewsFilter, TestTypes);

  TYPED_TEST(StdCxx20RangesViewsFilter, FilterEven)
  {
    TypeParam container = this->fill();
    int counter = 0;

    // Note: the simple forward implementation lacks the syntactic sugar of operator| to chain views
    for (const auto element : std_20::ranges::views::filter(
             std_20::ranges::views::all(container), [](int i) { return i % 2 == 0; }))
    {
      (void)element;
      counter++;
    }

    EXPECT_EQ(counter, 2);
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE
