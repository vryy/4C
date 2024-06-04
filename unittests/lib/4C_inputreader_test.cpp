/*---------------------------------------------------------------------------*/
/*! \file
\brief Unittests for reading input
\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_io_inputreader.hpp"

namespace
{
  using namespace FourC;

  TEST(ReadKeyValue, WithWhitespace)
  {
    const auto& [key, value] = CORE::IO::ReadKeyValue("key 1.0");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0");
  }

  TEST(ReadKeyValue, WithWhitespaceMultipleTakesFirst)
  {
    const auto& [key, value] = CORE::IO::ReadKeyValue("key 1.0 2.0 3");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0 2.0 3");
  }

  TEST(ReadKeyValue, WithWhitespaceAndEqualSignInside)
  {
    const auto& [key, value] = CORE::IO::ReadKeyValue("key=key value=value");
    EXPECT_EQ(key, "key=key");
    EXPECT_EQ(value, "value=value");
  }

  TEST(ReadKeyValue, WithEqualsSign)
  {
    const auto& [key, value] = CORE::IO::ReadKeyValue("key = 1.0");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0");
  }

  TEST(ReadKeyValue, WithEqualsSignMultipleTakesFirst)
  {
    const auto& [key, value] = CORE::IO::ReadKeyValue("key = 1.0 = 2.0 = 3.0=4.0");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0 = 2.0 = 3.0=4.0");
  }

  TEST(ReadKeyValue, WithEqualsSignNoKey) { EXPECT_ANY_THROW(CORE::IO::ReadKeyValue("   = 1.0")); }

  TEST(ReadKeyValue, WithEqualsSignNoKValue)
  {
    EXPECT_ANY_THROW(CORE::IO::ReadKeyValue(" key   =      "));
  }

  TEST(ReadKeyValue, SingleWordThrows) { EXPECT_ANY_THROW(CORE::IO::ReadKeyValue("key")); }


  TEST(ReadKeyValue, EmptyThrows) { EXPECT_ANY_THROW(CORE::IO::ReadKeyValue("")); };

  TEST(ReadKeyValue, WithEqualsSignNoSpaceThrows)
  {
    EXPECT_ANY_THROW(CORE::IO::ReadKeyValue("key=1.0"));
  }

}  // namespace