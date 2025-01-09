// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_utils_string.hpp"

#include <string>
#include <vector>

namespace
{
  using namespace FourC;
  using namespace Core::Utils;

  TEST(StringUtils, SplitStringSeparatorInString)
  {
    const std::string str = "1.1,2.2,3.3,4.4,5.5";
    auto split_string = split_string_list(str, ",");
    EXPECT_NO_THROW(split_string_list(str, ","));

    std::vector<std::string> expected_str = {"1.1", "2.2", "3.3", "4.4", "5.5"};
    EXPECT_EQ(split_string, expected_str);
  }

  TEST(StringUtils, StripComment)
  {
    const std::string& line = "abc//def";
    auto res = strip_comment(line);
    EXPECT_EQ(res, "abc");
  }

  TEST(StringUtils, StripCommentWithSpaceAfter)
  {
    const std::string& line = "abc// def";
    auto res = strip_comment(line);
    EXPECT_EQ(res, "abc");
  }

  TEST(StringUtils, StripCommentWithSpaceBefore)
  {
    const std::string& line = "abc //def";
    auto res = strip_comment(line);
    EXPECT_EQ(res, "abc");
  }
}  // namespace