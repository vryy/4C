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

}  // namespace