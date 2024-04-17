/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for LineDefinition

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_io_linedefinition.hpp"

#include <sstream>

namespace
{
  using namespace FourC;

  TEST(LineDefinitionTest, AddTag)
  {
    std::istringstream input("OMEGA");
    auto line_definition = INPUT::LineDefinition::Builder().AddTag("OMEGA").Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenTagRequiredButNothingGiven)
  {
    std::istringstream input("");
    auto line_definition = INPUT::LineDefinition::Builder().AddTag("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenTagRequiredButIntGiven)
  {
    std::istringstream input("1");
    auto line_definition = INPUT::LineDefinition::Builder().AddTag("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenTagRequiredButDoubleGiven)
  {
    std::istringstream input("1.23");
    auto line_definition = INPUT::LineDefinition::Builder().AddTag("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  // String
  TEST(LineDefinitionTest, AddString)
  {
    std::istringstream input("OMEGA");
    auto line_definition = INPUT::LineDefinition::Builder().AddString("OMEGA").Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenStringRequiredButNothingGiven)
  {
    std::istringstream input("");
    auto line_definition = INPUT::LineDefinition::Builder().AddString("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  // Int
  TEST(LineDefinitionTest, AddInt)
  {
    std::istringstream input("1");
    auto line_definition = INPUT::LineDefinition::Builder().AddInt("OMEGA").Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntRequiredButNothingGiven)
  {
    std::istringstream input("");
    auto line_definition = INPUT::LineDefinition::Builder().AddInt("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntRequiredButDoubleGiven)
  {
    std::istringstream input("1.23");
    auto line_definition = INPUT::LineDefinition::Builder().AddInt("OMEGA").Build();
    EXPECT_ANY_THROW(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntConcatenatedWithSomething)
  {
    std::istringstream input("1*8e+2");
    auto line_definition = INPUT::LineDefinition::Builder().AddInt("OMEGA").Build();
    EXPECT_ANY_THROW(line_definition.Read(input));
  }

  // Int Vector
  TEST(LineDefinitionTest, AddIntVector)
  {
    std::istringstream input("1 2 3");
    auto line_definition = INPUT::LineDefinition::Builder().AddIntVector("OMEGA", 3).Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("1 2");
    auto line_definition = INPUT::LineDefinition::Builder().AddIntVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("1 2 3 4");
    auto line_definition = INPUT::LineDefinition::Builder().AddIntVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenIntVectorRequiredButDoubleVectorEntriesGiven)
  {
    std::istringstream input("1.23 2.34 3.45");
    auto line_definition = INPUT::LineDefinition::Builder().AddIntVector("OMEGA", 3).Build();
    EXPECT_ANY_THROW(line_definition.Read(input));
  }

  // Double Vector
  TEST(LineDefinitionTest, AddDoubleVector)
  {
    std::istringstream input("1 2 3");
    auto line_definition = INPUT::LineDefinition::Builder().AddDoubleVector("OMEGA", 3).Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenDoubleVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("1 2");
    auto line_definition = INPUT::LineDefinition::Builder().AddDoubleVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenDoubleVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("1 2 3 4");
    auto line_definition = INPUT::LineDefinition::Builder().AddDoubleVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  // Named String
  TEST(LineDefinitionTest, AddNamedString)
  {
    std::istringstream input("OMEGA TEST");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedString("OMEGA").Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedStringRequiredButNothingGiven)
  {
    std::istringstream input("OMEGA ");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedString("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedStringRequiredButNoNameGiven)
  {
    std::istringstream input("TEST");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedString("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  // Named Int
  TEST(LineDefinitionTest, AddNamedInt)
  {
    std::istringstream input("OMEGA 1");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedInt("OMEGA").Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntRequiredButNothingGiven)
  {
    std::istringstream input("OMEGA");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedInt("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntRequiredButDoubleGiven)
  {
    std::istringstream input("OMEGA 1.23");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedInt("OMEGA").Build();
    EXPECT_ANY_THROW(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntRequiredButNoNameGiven)
  {
    std::istringstream input("1");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedInt("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  // Named Int Vector
  TEST(LineDefinitionTest, AddNamedIntVector)
  {
    std::istringstream input("OMEGA 1 2 3");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedIntVector("OMEGA", 3).Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1 2");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedIntVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1 2 3 4");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedIntVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButDoubleVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedIntVector("OMEGA", 3).Build();
    EXPECT_ANY_THROW(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButNoNameGiven)
  {
    std::istringstream input("1.23 2.34 3.45");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedIntVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  // Named Double
  TEST(LineDefinitionTest, AddNamedDouble)
  {
    std::istringstream input("OMEGA 1.23");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedDouble("OMEGA").Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleRequiredButNothingGiven)
  {
    std::istringstream input("OMEGA");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedDouble("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleConcatenatedWithSomething)
  {
    std::istringstream input("OMEGA 123.45*893");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedDouble("OMEGA").Build();
    EXPECT_ANY_THROW(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleRequiredButNoNameGiven)
  {
    std::istringstream input("1.23");
    auto line_definition = INPUT::LineDefinition::Builder().AddNamedDouble("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  // Named Double Vector
  TEST(LineDefinitionTest, AddNamedDoubleVector)
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    auto line_definition =
        INPUT::LineDefinition::Builder().AddNamedDoubleVector("OMEGA", 3).Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1.23 2.34");
    auto line_definition =
        INPUT::LineDefinition::Builder().AddNamedDoubleVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45 4.56");
    auto line_definition =
        INPUT::LineDefinition::Builder().AddNamedDoubleVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadFalseWhenNamedDoubleVectorRequiredButNoNameGiven)
  {
    std::istringstream input("1.23 2.34 3.45");
    auto line_definition =
        INPUT::LineDefinition::Builder().AddNamedDoubleVector("OMEGA", 3).Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  // Too much in stream except from comments and whitespaces
  TEST(LineDefinitionTest, ReadFalseWhenStreamHasTooMuch)
  {
    std::istringstream input("OMEGA More data that is not in lineDefinition.");
    auto line_definition = INPUT::LineDefinition::Builder().AddTag("OMEGA").Build();
    EXPECT_FALSE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadTrueWhenStreamHasTooMuchButOnlyComments)
  {
    std::istringstream input("OMEGA // Comment");
    auto line_definition = INPUT::LineDefinition::Builder().AddTag("OMEGA").Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionTest, ReadTrue_When_StreamHasTooMuchButOnlyWhitespaces)
  {
    std::istringstream input("OMEGA           ");
    auto line_definition = INPUT::LineDefinition::Builder().AddTag("OMEGA").Build();
    EXPECT_TRUE(line_definition.Read(input));
  }

  // Empty stream
  TEST(LineDefinitionTest, ReadTrueWhenEmptyStreamIntoEmptyDefinition)
  {
    std::istringstream input("");
    INPUT::LineDefinition line_definition;
    EXPECT_TRUE(line_definition.Read(input));
  }

  TEST(LineDefinitionPrinting, Empty)
  {
    std::ostringstream out;
    INPUT::LineDefinition().Print(out);
    EXPECT_EQ(out.str(), "");
  }

  TEST(LineDefinitionPrinting, ComplicatedLine)
  {
    std::ostringstream out;
    INPUT::LineDefinition::Builder()
        .AddTag("abc")
        .AddInt("i")
        .AddNamedDouble("d")
        .AddNamedIntVector("iv", 3)
        .AddOptionalNamedStringVector("s", 2)
        .AddOptionalNamedPairOfStringAndDoubleVector("pairs", INPUT::LengthFromIntNamed("i"))
        .Build()
        .Print(out);

    EXPECT_EQ(out.str(), "abc 0 d 0 iv 0 0 0  [ pairs [...] s '' ''  ] ");
  }
}  // namespace
