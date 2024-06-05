/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for functions in 4C_io_utils_reader.H

\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_utils_string.hpp"

#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
  using namespace FourC;
  using namespace CORE::UTILS;

  TEST(UtilsReaderTests, SplitStringSeparatorInString)
  {
    const std::string str = "1.1,2.2,3.3,4.4,5.5";
    auto splitted_string = SplitStringList(str, ",");
    EXPECT_NO_THROW(SplitStringList(str, ","));

    std::vector<std::string> expected_str = {"1.1", "2.2", "3.3", "4.4", "5.5"};
    EXPECT_EQ(splitted_string, expected_str);
  }
}  // namespace