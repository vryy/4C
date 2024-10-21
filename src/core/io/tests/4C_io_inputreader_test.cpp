// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_io_inputreader.hpp"

namespace
{
  using namespace FourC;

  TEST(ReadKeyValue, WithWhitespace)
  {
    const auto& [key, value] = Core::IO::read_key_value("key 1.0");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0");
  }

  TEST(ReadKeyValue, WithWhitespaceMultipleTakesFirst)
  {
    const auto& [key, value] = Core::IO::read_key_value("key 1.0 2.0 3");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0 2.0 3");
  }

  TEST(ReadKeyValue, WithWhitespaceAndEqualSignInside)
  {
    const auto& [key, value] = Core::IO::read_key_value("key=key value=value");
    EXPECT_EQ(key, "key=key");
    EXPECT_EQ(value, "value=value");
  }

  TEST(ReadKeyValue, WithEqualsSign)
  {
    const auto& [key, value] = Core::IO::read_key_value("key = 1.0");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0");
  }

  TEST(ReadKeyValue, WithEqualsSignMultipleTakesFirst)
  {
    const auto& [key, value] = Core::IO::read_key_value("key = 1.0 = 2.0 = 3.0=4.0");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0 = 2.0 = 3.0=4.0");
  }

  TEST(ReadKeyValue, WithEqualsSignNoKey)
  {
    EXPECT_ANY_THROW(Core::IO::read_key_value("   = 1.0"));
  }

  TEST(ReadKeyValue, WithEqualsSignNoKValue)
  {
    EXPECT_ANY_THROW(Core::IO::read_key_value(" key   =      "));
  }

  TEST(ReadKeyValue, SingleWordThrows) { EXPECT_ANY_THROW(Core::IO::read_key_value("key")); }


  TEST(ReadKeyValue, EmptyThrows) { EXPECT_ANY_THROW(Core::IO::read_key_value("")); };

  TEST(ReadKeyValue, WithEqualsSignNoSpaceThrows)
  {
    EXPECT_ANY_THROW(Core::IO::read_key_value("key=1.0"));
  }

  TEST(StreamLineIterator, Empty)
  {
    auto stream = std::make_shared<std::istringstream>("");
    Core::IO::Internal::StreamLineIterator it{stream};
    Core::IO::Internal::StreamLineIterator it_end{};
    EXPECT_EQ(it, it_end);
  }

  TEST(StreamLineIterator, SingleLine)
  {
    auto stream = std::make_shared<std::istringstream>("test");
    Core::IO::Internal::StreamLineIterator it{stream};
    Core::IO::Internal::StreamLineIterator it_end{};
    std::string line;
    for (; it != it_end; ++it)
    {
      line += *it;
    }
    EXPECT_EQ(line, "test");
  }

  TEST(StreamLineIterator, MultipleLineUntilEnd)
  {
    auto stream = std::make_shared<std::istringstream>("a\nb\nc\n");
    Core::IO::Internal::StreamLineIterator it{stream};
    Core::IO::Internal::StreamLineIterator it_end{};
    std::string line;
    for (; it != it_end; ++it)
    {
      line += *it;
    }
    EXPECT_EQ(line, "abc");
  }

  TEST(StreamLineIterator, MultipleLineUntilGivenLine)
  {
    auto stream = std::make_shared<std::istringstream>("a\nb\nc\n");
    Core::IO::Internal::StreamLineIterator it{stream, 2};
    Core::IO::Internal::StreamLineIterator it_end{};
    std::string line;
    for (; it != it_end; ++it)
    {
      line += *it;
    }
    EXPECT_EQ(line, "ab");
  }

  TEST(StreamLineIterator, EmptyRange)
  {
    auto stream = std::make_shared<std::istringstream>("a\nb\nc\n");
    Core::IO::Internal::StreamLineIterator it{stream, 0};
    // Read zero lines
    Core::IO::Internal::StreamLineIterator it_end{};
    std::string line;
    for (; it != it_end; ++it)
    {
      line += *it;
    }
    EXPECT_EQ(line, "");
  }

}  // namespace