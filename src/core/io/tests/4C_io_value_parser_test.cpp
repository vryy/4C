/*---------------------------------------------------------------------------*/
/*! \file
\brief Unittests for value parser
\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_io_value_parser.hpp"
namespace
{

  using namespace FourC;

  TEST(ValueParser, ConsumeSuccess)
  {
    std::istringstream string_stream("expected");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");
    parser.consume("expected");
  }

  TEST(ValueParser, ConsumeFail)
  {
    std::istringstream string_stream("unexpected");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");
    EXPECT_ANY_THROW(parser.consume("expected"));
  }

  TEST(ValueParser, ReadDoubleSuccess)
  {
    std::istringstream string_stream("11.3");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    EXPECT_EQ(parser.read<double>(), 11.3);
  }

  TEST(ValueParser, ReadIntSuccess)
  {
    std::istringstream string_stream("42");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    EXPECT_EQ(parser.read<int>(), 42);
  }

  TEST(ValueParser, ReadDoubleFromIntSuccess)
  {
    std::istringstream string_stream("42");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    EXPECT_EQ(parser.read<double>(), 42);
  }

  TEST(ValueParser, ReadIntFromDoubleFail)
  {
    std::istringstream string_stream("11.3");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    EXPECT_ANY_THROW(parser.read<int>());
  }

  TEST(ValueParser, ReadFailWithExtraCharacters)
  {
    std::istringstream string_stream("11.3*2");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    EXPECT_ANY_THROW(parser.read<double>());
  }

  TEST(ValueParser, ReadExtraCharactersAfterWhitespaceCheckEOF)
  {
    std::istringstream string_stream("11.3 3");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    EXPECT_EQ(parser.read<double>(), 11.3);
    EXPECT_FALSE(parser.eof());
  }

  TEST(ValueParser, ReadIntArraySuccess)
  {
    std::istringstream string_stream("1 2 3");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    const auto array = parser.read_array<int, 3>();

    EXPECT_EQ(array[0], 1);
    EXPECT_EQ(array[1], 2);
    EXPECT_EQ(array[2], 3);
    EXPECT_TRUE(parser.eof());
  }

  TEST(ValueParser, ReadIntArrayFailTooLong)
  {
    std::istringstream string_stream("1 2 3 4");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    parser.read_array<int, 3>();
    EXPECT_FALSE(parser.eof());
  }

  TEST(ValueParser, ReadIntArrayFailTooShort)
  {
    std::istringstream string_stream("1 2");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    EXPECT_ANY_THROW((parser.read_array<int, 3>()));
  }

  TEST(ValueParser, ReadIntArrayFailOtherCharacters)
  {
    std::istringstream string_stream("1 2 a");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    EXPECT_ANY_THROW((parser.read_array<int, 3>()));
  }

  TEST(ValueParser, ReadCombinedIntStringArraySuccess)
  {
    std::istringstream string_stream("1 2 3 a b c");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    const auto ints = parser.read_array<int, 3>();
    const auto strings = parser.read_array<std::string, 3>();

    EXPECT_EQ(ints[0], 1);
    EXPECT_EQ(ints[1], 2);
    EXPECT_EQ(ints[2], 3);
    EXPECT_EQ(strings[0], "a");
    EXPECT_EQ(strings[1], "b");
    EXPECT_EQ(strings[2], "c");
    EXPECT_TRUE(parser.eof());
  }

  TEST(ValueParser, ReadDoubleVectorSuccess)
  {
    // read a vector of doubles with not specified size
    std::istringstream string_stream("0.1 0.2 0.3 0.4 0.5 0.6");
    Core::IO::ValueParser parser(string_stream, "While reading section MY PARAMETERS: ");

    std::vector<double> vec;

    while (!parser.eof())
    {
      vec.push_back(parser.read<double>());
    }

    for (size_t i = 0; i < 6; ++i)
    {
      EXPECT_DOUBLE_EQ(vec[i], 0.1 * (i + 1));
    }
  }


}  // namespace