/*---------------------------------------------------------------------------*/
/*! \file
\brief Unittests for reading input
\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_io_inputreader.H"

namespace
{
  using namespace BACI;

  TEST(ReadKeyValue, WithWhitespace)
  {
    const auto& [key, value] = DRT::INPUT::ReadKeyValue("key 1.0");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0");
  }

  TEST(ReadKeyValue, WithWhitespaceMultipleTakesFirst)
  {
    const auto& [key, value] = DRT::INPUT::ReadKeyValue("key 1.0 2.0 3");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0 2.0 3");
  }

  TEST(ReadKeyValue, WithWhitespaceAndEqualSignInside)
  {
    const auto& [key, value] = DRT::INPUT::ReadKeyValue("key=key value=value");
    EXPECT_EQ(key, "key=key");
    EXPECT_EQ(value, "value=value");
  }

  TEST(ReadKeyValue, WithEqualsSign)
  {
    const auto& [key, value] = DRT::INPUT::ReadKeyValue("key = 1.0");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0");
  }

  TEST(ReadKeyValue, WithEqualsSignMultipleTakesFirst)
  {
    const auto& [key, value] = DRT::INPUT::ReadKeyValue("key = 1.0 = 2.0 = 3.0=4.0");
    EXPECT_EQ(key, "key");
    EXPECT_EQ(value, "1.0 = 2.0 = 3.0=4.0");
  }

  TEST(ReadKeyValue, WithEqualsSignNoKey)
  {
    EXPECT_ANY_THROW(DRT::INPUT::ReadKeyValue("   = 1.0"));
  }

  TEST(ReadKeyValue, WithEqualsSignNoKValue)
  {
    EXPECT_ANY_THROW(DRT::INPUT::ReadKeyValue(" key   =      "));
  }

  TEST(ReadKeyValue, SingleWordThrows) { EXPECT_ANY_THROW(DRT::INPUT::ReadKeyValue("key")); }


  TEST(ReadKeyValue, EmptyThrows) { EXPECT_ANY_THROW(DRT::INPUT::ReadKeyValue("")); };

  TEST(ReadKeyValue, WithEqualsSignNoSpaceThrows)
  {
    EXPECT_ANY_THROW(DRT::INPUT::ReadKeyValue("key=1.0"));
  }

}  // namespace