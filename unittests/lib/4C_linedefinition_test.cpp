/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for LineDefinition

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_io_linedefinition.hpp"

#include <sstream>

namespace
{
  using namespace FourC;

  TEST(LineDefinitionTest, add_tag)
  {
    std::istringstream input("OMEGA");
    auto line_definition = Input::LineDefinition::Builder().add_tag("OMEGA").build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenTagRequiredButNothingGiven)
  {
    std::istringstream input("");
    auto line_definition = Input::LineDefinition::Builder().add_tag("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenTagRequiredButIntGiven)
  {
    std::istringstream input("1");
    auto line_definition = Input::LineDefinition::Builder().add_tag("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenTagRequiredButDoubleGiven)
  {
    std::istringstream input("1.23");
    auto line_definition = Input::LineDefinition::Builder().add_tag("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  // String
  TEST(LineDefinitionTest, add_string)
  {
    std::istringstream input("OMEGA");
    auto line_definition = Input::LineDefinition::Builder().add_string("OMEGA").build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenStringRequiredButNothingGiven)
  {
    std::istringstream input("");
    auto line_definition = Input::LineDefinition::Builder().add_string("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  // Int
  TEST(LineDefinitionTest, add_int)
  {
    std::istringstream input("1");
    auto line_definition = Input::LineDefinition::Builder().add_int("OMEGA").build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntRequiredButNothingGiven)
  {
    std::istringstream input("");
    auto line_definition = Input::LineDefinition::Builder().add_int("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntRequiredButDoubleGiven)
  {
    std::istringstream input("1.23");
    auto line_definition = Input::LineDefinition::Builder().add_int("OMEGA").build();
    EXPECT_ANY_THROW(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntConcatenatedWithSomething)
  {
    std::istringstream input("1*8e+2");
    auto line_definition = Input::LineDefinition::Builder().add_int("OMEGA").build();
    EXPECT_ANY_THROW(line_definition.read(input));
  }

  // Int Vector
  TEST(LineDefinitionTest, add_int_vector)
  {
    std::istringstream input("1 2 3");
    auto line_definition = Input::LineDefinition::Builder().add_int_vector("OMEGA", 3).build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("1 2");
    auto line_definition = Input::LineDefinition::Builder().add_int_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("1 2 3 4");
    auto line_definition = Input::LineDefinition::Builder().add_int_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntVectorRequiredButDoubleVectorEntriesGiven)
  {
    std::istringstream input("1.23 2.34 3.45");
    auto line_definition = Input::LineDefinition::Builder().add_int_vector("OMEGA", 3).build();
    EXPECT_ANY_THROW(line_definition.read(input));
  }

  // Double Vector
  TEST(LineDefinitionTest, add_double_vector)
  {
    std::istringstream input("1 2 3");
    auto line_definition = Input::LineDefinition::Builder().add_double_vector("OMEGA", 3).build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenDoubleVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("1 2");
    auto line_definition = Input::LineDefinition::Builder().add_double_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenDoubleVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("1 2 3 4");
    auto line_definition = Input::LineDefinition::Builder().add_double_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  // Named String
  TEST(LineDefinitionTest, add_named_string)
  {
    std::istringstream input("OMEGA TEST");
    auto line_definition = Input::LineDefinition::Builder().add_named_string("OMEGA").build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedStringRequiredButNothingGiven)
  {
    std::istringstream input("OMEGA ");
    auto line_definition = Input::LineDefinition::Builder().add_named_string("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedStringRequiredButNoNameGiven)
  {
    std::istringstream input("TEST");
    auto line_definition = Input::LineDefinition::Builder().add_named_string("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  // Named Int
  TEST(LineDefinitionTest, add_named_int)
  {
    std::istringstream input("OMEGA 1");
    auto line_definition = Input::LineDefinition::Builder().add_named_int("OMEGA").build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntRequiredButNothingGiven)
  {
    std::istringstream input("OMEGA");
    auto line_definition = Input::LineDefinition::Builder().add_named_int("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntRequiredButDoubleGiven)
  {
    std::istringstream input("OMEGA 1.23");
    auto line_definition = Input::LineDefinition::Builder().add_named_int("OMEGA").build();
    EXPECT_ANY_THROW(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntRequiredButNoNameGiven)
  {
    std::istringstream input("1");
    auto line_definition = Input::LineDefinition::Builder().add_named_int("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  // Named Int Vector
  TEST(LineDefinitionTest, add_named_int_vector)
  {
    std::istringstream input("OMEGA 1 2 3");
    auto line_definition =
        Input::LineDefinition::Builder().add_named_int_vector("OMEGA", 3).build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1 2");
    auto line_definition =
        Input::LineDefinition::Builder().add_named_int_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1 2 3 4");
    auto line_definition =
        Input::LineDefinition::Builder().add_named_int_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButDoubleVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    auto line_definition =
        Input::LineDefinition::Builder().add_named_int_vector("OMEGA", 3).build();
    EXPECT_ANY_THROW(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButNoNameGiven)
  {
    std::istringstream input("1.23 2.34 3.45");
    auto line_definition =
        Input::LineDefinition::Builder().add_named_int_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  // Named Double
  TEST(LineDefinitionTest, add_named_double)
  {
    std::istringstream input("OMEGA 1.23");
    auto line_definition = Input::LineDefinition::Builder().add_named_double("OMEGA").build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleRequiredButNothingGiven)
  {
    std::istringstream input("OMEGA");
    auto line_definition = Input::LineDefinition::Builder().add_named_double("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleConcatenatedWithSomething)
  {
    std::istringstream input("OMEGA 123.45*893");
    auto line_definition = Input::LineDefinition::Builder().add_named_double("OMEGA").build();
    EXPECT_ANY_THROW(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleRequiredButNoNameGiven)
  {
    std::istringstream input("1.23");
    auto line_definition = Input::LineDefinition::Builder().add_named_double("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  // Named Double Vector
  TEST(LineDefinitionTest, add_named_double_vector)
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    auto line_definition =
        Input::LineDefinition::Builder().add_named_double_vector("OMEGA", 3).build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1.23 2.34");
    auto line_definition =
        Input::LineDefinition::Builder().add_named_double_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45 4.56");
    auto line_definition =
        Input::LineDefinition::Builder().add_named_double_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleVectorRequiredButNoNameGiven)
  {
    std::istringstream input("1.23 2.34 3.45");
    auto line_definition =
        Input::LineDefinition::Builder().add_named_double_vector("OMEGA", 3).build();
    EXPECT_FALSE(line_definition.read(input));
  }

  // Too much in stream except from comments and whitespaces
  TEST(LineDefinitionTest, ReadFalseWhenStreamHasTooMuch)
  {
    std::istringstream input("OMEGA More data that is not in lineDefinition.");
    auto line_definition = Input::LineDefinition::Builder().add_tag("OMEGA").build();
    EXPECT_FALSE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadTrueWhenStreamHasTooMuchButOnlyComments)
  {
    std::istringstream input("OMEGA // Comment");
    auto line_definition = Input::LineDefinition::Builder().add_tag("OMEGA").build();
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionTest, ReadTrue_When_StreamHasTooMuchButOnlyWhitespaces)
  {
    std::istringstream input("OMEGA           ");
    auto line_definition = Input::LineDefinition::Builder().add_tag("OMEGA").build();
    EXPECT_TRUE(line_definition.read(input));
  }

  // Empty stream
  TEST(LineDefinitionTest, ReadTrueWhenEmptyStreamIntoEmptyDefinition)
  {
    std::istringstream input("");
    Input::LineDefinition line_definition;
    EXPECT_TRUE(line_definition.read(input));
  }

  TEST(LineDefinitionPrinting, Empty)
  {
    std::ostringstream out;
    Input::LineDefinition().print(out);
    EXPECT_EQ(out.str(), "");
  }

  TEST(LineDefinitionPrinting, ComplicatedLine)
  {
    std::ostringstream out;
    Input::LineDefinition::Builder()
        .add_tag("abc")
        .add_int("i")
        .add_named_double("d")
        .add_named_int_vector("iv", 3)
        .add_optional_named_string_vector("s", 2)
        .add_optional_named_pair_of_string_and_double_vector(
            "pairs", Input::LengthFromIntNamed("i"))
        .build()
        .print(out);

    EXPECT_EQ(out.str(), "abc 0 d 0 iv 0 0 0  [ pairs [...] s '' ''  ] ");
  }
}  // namespace
