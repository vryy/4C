/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for the forward implementation of the range functionality in
utils_std_cxx20_ranges.H

\level 0
*/

#include <gtest/gtest.h>

#include "baci_utils_std_cxx20_ranges.H"

#include <list>

namespace
{
  template <typename Container>
  class StdCxx20RangesViewsAll : public testing::Test
  {
   protected:
    [[nodiscard]] Container Fill() const { return Container{1, 2, 3, 4, 5}; }
  };


  using TestTypes = ::testing::Types<std::array<int, 5>, std::vector<int>, std::list<int>>;
  TYPED_TEST_SUITE(StdCxx20RangesViewsAll, TestTypes);

  TYPED_TEST(StdCxx20RangesViewsAll, testIterator)
  {
    TypeParam container = this->Fill();
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
    const TypeParam container = this->Fill();
    int sum_of_elements = 0;

    for (const auto element : std_20::ranges::views::all(container))
    {
      sum_of_elements += element;
    }

    EXPECT_EQ(sum_of_elements, 15);
  }

  TYPED_TEST(StdCxx20RangesViewsAll, testEleManipulation)
  {
    TypeParam container = this->Fill();

    for (auto &element : std_20::ranges::views::all(container))
    {
      element += element;
    }

    EXPECT_EQ(container, (TypeParam{2, 4, 6, 8, 10}));
  }
}  // namespace