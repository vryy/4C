// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include <gmock/gmock.h>

#include "4C_utils_vector2D.hpp"

namespace
{
  using namespace FourC::Core::Utils;
  using ::testing::ElementsAre;

  template <typename T, std::size_t cols>
  Vector2D<T> fill_vector(std::vector<std::array<T, cols>> data)
  {
    Vector2D<T> v{cols};
    v.reserve(data.size());

    for (const auto& row : data)
    {
      std::vector<T> v_row(row.begin(), row.end());
      v.push_back(v_row);
    }
    return v;
  }


  TEST(Vector2DTest, Default)
  {
    const Vector2D<int> v;
    EXPECT_EQ(v.size(), 0);
    EXPECT_EQ(v.num_components(), 0);
  }

  TEST(Vector2DTest, ReserveAndPushBack)
  {
    Vector2D<int> v = fill_vector<int, 3>({{1, 2, 3}, {4, 5, 6}});

    ASSERT_THAT(v.at(0), ElementsAre(1, 2, 3));
    ASSERT_THAT(v.at(1), ElementsAre(4, 5, 6));
  }

  TEST(Vector2DTest, AccessOperatorRead)
  {
    Vector2D<int> v = fill_vector<int, 3>({{1, 2, 3}, {4, 5, 6}});

    EXPECT_EQ(v(0, 0), 1);
    EXPECT_EQ(v(0, 1), 2);
    EXPECT_EQ(v(0, 2), 3);
    EXPECT_EQ(v(1, 0), 4);
    EXPECT_EQ(v(1, 1), 5);
    EXPECT_EQ(v(1, 2), 6);
  }

  TEST(Vector2DTest, AccessOperatorWrite)
  {
    Vector2D<int> v = fill_vector<int, 3>({{1, 2, 3}, {4, 5, 6}});
    v(0, 1) = 3;

    EXPECT_EQ(v(0, 0), 1);
    EXPECT_EQ(v(0, 1), 3);
    EXPECT_EQ(v(0, 2), 3);
    EXPECT_EQ(v(1, 0), 4);
    EXPECT_EQ(v(1, 1), 5);
    EXPECT_EQ(v(1, 2), 6);
  }

  TEST(Vector2DTest, AccessRow)
  {
    Vector2D<bool> v = fill_vector<bool, 3>({{true, false, true}, {false, true, false}});

    ASSERT_THAT(v.at(0), ElementsAre(true, false, true));
    ASSERT_THAT(v.at(1), ElementsAre(false, true, false));
  }

  TEST(Vector2DTest, Rows)
  {
    Vector2D<bool> v = fill_vector<bool, 3>({{true, false, true}, {false, true, false}});

    std::size_t i = 0;
    for (const auto row : v.items())
    {
      if (i == 0)
        ASSERT_THAT(row, ElementsAre(true, false, true));
      else if (i == 1)
        ASSERT_THAT(row, ElementsAre(false, true, false));
      else
        FAIL() << "Unexpected row index " << i;
      ++i;
    }
  }

}  // namespace
