// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_utils_lazy_ptr.hpp"

namespace
{
  using namespace FourC::Core::Utils;

  struct MoveOnlyType
  {
    MoveOnlyType() = default;
    MoveOnlyType(const MoveOnlyType&) = delete;
    MoveOnlyType(MoveOnlyType&&) = default;
    MoveOnlyType& operator=(const MoveOnlyType&) = delete;
    MoveOnlyType& operator=(MoveOnlyType&&) = default;

    int value = 42;
  };

  TEST(LazyPtr, BasicConstruction)
  {
    const LazyPtr<int> a([]() { return std::make_unique<int>(42); });
    const LazyPtr<int> b([]() { return std::make_unique<int>(42); });
    const LazyPtr<int> c = a;
    EXPECT_EQ(*a + *b + *c, 126);
  }

  TEST(LazyPtr, ConstructFromSharedUnique)
  {
    const LazyPtr<int> a(std::make_shared<int>(42));
    const LazyPtr<int> b(std::make_unique<int>(42));
    const LazyPtr<int> c = a;
    EXPECT_EQ(*a + *b + *c, 126);
  }

  TEST(LazyPtr, MoveOnlyType)
  {
    const LazyPtr<MoveOnlyType> a([]() { return std::make_unique<MoveOnlyType>(); });
    const LazyPtr<MoveOnlyType> b([]() { return std::make_unique<MoveOnlyType>(); });
    EXPECT_EQ(a->value + b->value, 84);
  }

  TEST(LazyPtr, DependencyChain)
  {
    std::map<int, LazyPtr<int>> ptr_list{
        {0, LazyPtr<int>([&ptr_list]() { return std::make_shared<int>(*ptr_list.at(1) + 1); })},
        {1, LazyPtr<int>([]() { return std::make_unique<int>(2); })},
    };

    EXPECT_EQ(*ptr_list.at(0), 3);
  }
}  // namespace
