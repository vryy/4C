/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for LineDefinition

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_lib_linedefinition.H"

namespace
{
  class LineDefinitionTest : public testing::Test
  {
   public:
    LineDefinitionTest() { lineDefinition = Teuchos::rcp(new DRT::INPUT::LineDefinition()); }

    void ExpectReadReturnsTrue(std::istringstream& input)
    {
      EXPECT_TRUE(lineDefinition->Read(input));
    }

    void ExpectReadReturnsFalse(std::istringstream& input)
    {
      EXPECT_TRUE(not lineDefinition->Read(input));
    }

    void ExpectThrowsError(std::istringstream& input)
    {
      EXPECT_ANY_THROW(lineDefinition->Read(input));
    }

   protected:
    Teuchos::RCP<DRT::INPUT::LineDefinition> lineDefinition;
  };

  TEST_F(LineDefinitionTest, AddTag)
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddTag("OMEGA");
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenTagRequiredButNothingGiven)
  {
    std::istringstream input("");
    lineDefinition->AddTag("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenTagRequiredButIntGiven)
  {
    std::istringstream input("1");
    lineDefinition->AddTag("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenTagRequiredButDoubleGiven)
  {
    std::istringstream input("1.23");
    lineDefinition->AddTag("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  // String
  TEST_F(LineDefinitionTest, AddString)
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddString("OMEGA");
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenStringRequiredButNothingGiven)
  {
    std::istringstream input("");
    lineDefinition->AddString("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  // Int
  TEST_F(LineDefinitionTest, AddInt)
  {
    std::istringstream input("1");
    lineDefinition->AddInt("OMEGA");
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenIntRequiredButNothingGiven)
  {
    std::istringstream input("");
    lineDefinition->AddInt("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenIntRequiredButDoubleGiven)
  {
    std::istringstream input("1.23");
    lineDefinition->AddInt("OMEGA");
    ExpectThrowsError(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenIntConcatenatedWithSomething)
  {
    std::istringstream input("1*8e+2");
    lineDefinition->AddInt("OMEGA");
    ExpectThrowsError(input);
  }

  // Int Vector
  TEST_F(LineDefinitionTest, AddIntVector)
  {
    std::istringstream input("1 2 3");
    lineDefinition->AddIntVector("OMEGA", 3);
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenIntVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("1 2");
    lineDefinition->AddIntVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenIntVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("1 2 3 4");
    lineDefinition->AddIntVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenIntVectorRequiredButDoubleVectorEntriesGiven)
  {
    std::istringstream input("1.23 2.34 3.45");
    lineDefinition->AddIntVector("OMEGA", 3);
    ExpectThrowsError(input);
  }

  // Double Vector
  TEST_F(LineDefinitionTest, AddDoubleVector)
  {
    std::istringstream input("1 2 3");
    lineDefinition->AddDoubleVector("OMEGA", 3);
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenDoubleVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("1 2");
    lineDefinition->AddDoubleVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenDoubleVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("1 2 3 4");
    lineDefinition->AddDoubleVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  // Named String
  TEST_F(LineDefinitionTest, AddNamedString)
  {
    std::istringstream input("OMEGA TEST");
    lineDefinition->AddNamedString("OMEGA");
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedStringRequiredButNothingGiven)
  {
    std::istringstream input("OMEGA ");
    lineDefinition->AddNamedString("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedStringRequiredButNoNameGiven)
  {
    std::istringstream input("TEST");
    lineDefinition->AddNamedString("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  // Named Int
  TEST_F(LineDefinitionTest, AddNamedInt)
  {
    std::istringstream input("OMEGA 1");
    lineDefinition->AddNamedInt("OMEGA");
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedIntRequiredButNothingGiven)
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddNamedInt("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedIntRequiredButDoubleGiven)
  {
    std::istringstream input("OMEGA 1.23");
    lineDefinition->AddNamedInt("OMEGA");
    ExpectThrowsError(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedIntRequiredButNoNameGiven)
  {
    std::istringstream input("1");
    lineDefinition->AddNamedInt("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  // Named Int Vector
  TEST_F(LineDefinitionTest, AddNamedIntVector)
  {
    std::istringstream input("OMEGA 1 2 3");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1 2");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1 2 3 4");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButDoubleVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    ExpectThrowsError(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedIntVectorRequiredButNoNameGiven)
  {
    std::istringstream input("1.23 2.34 3.45");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  // Named Double
  TEST_F(LineDefinitionTest, AddNamedDouble)
  {
    std::istringstream input("OMEGA 1.23");
    lineDefinition->AddNamedDouble("OMEGA");
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedDoubleRequiredButNothingGiven)
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddNamedDouble("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedDoubleConcatenatedWithSomething)
  {
    std::istringstream input("OMEGA 123.45*893");
    lineDefinition->AddNamedDouble("OMEGA");
    ExpectThrowsError(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedDoubleRequiredButNoNameGiven)
  {
    std::istringstream input("1.23");
    lineDefinition->AddNamedDouble("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  // Named Double Vector
  TEST_F(LineDefinitionTest, AddNamedDoubleVector)
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedDoubleVectorRequiredButTooFewVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1.23 2.34");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedDoubleVectorRequiredButTooManyVectorEntriesGiven)
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45 4.56");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadFalseWhenNamedDoubleVectorRequiredButNoNameGiven)
  {
    std::istringstream input("1.23 2.34 3.45");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    ExpectReadReturnsFalse(input);
  }

  // Too much in stream except from comments and whitespaces
  TEST_F(LineDefinitionTest, ReadFalseWhenStreamHasTooMuch)
  {
    std::istringstream input("OMEGA More data that is not in lineDefinition.");
    lineDefinition->AddTag("OMEGA");
    ExpectReadReturnsFalse(input);
  }

  TEST_F(LineDefinitionTest, ReadTrueWhenStreamHasTooMuchButOnlyComments)
  {
    std::istringstream input("OMEGA // Comment");
    lineDefinition->AddTag("OMEGA");
    ExpectReadReturnsTrue(input);
  }

  TEST_F(LineDefinitionTest, ReadTrue_When_StreamHasTooMuchButOnlyWhitespaces)
  {
    std::istringstream input("OMEGA           ");
    lineDefinition->AddTag("OMEGA");
    ExpectReadReturnsTrue(input);
  }

  // Empty stream
  TEST_F(LineDefinitionTest, ReadTrueWhenEmptyStreamIntoEmptyDefinition)
  {
    std::istringstream input("");
    ExpectReadReturnsTrue(input);
  }
}  // namespace