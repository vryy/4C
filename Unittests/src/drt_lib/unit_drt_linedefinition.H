/*! \file

\brief unit tests for LineDefinition

\level 3
*/
#ifndef UNIT_DRT_LINEDEFINITION_H_
#define UNIT_DRT_LINEDEFINITION_H_

#include <cxxtest/TestSuite.h>
#include "src/drt_lib/drt_linedefinition.H"

namespace DRT
{
  namespace INPUT
  {
    class LineDefinition_TestSuite;
  }
}  // namespace DRT


class DRT::INPUT::LineDefinition_TestSuite : public CxxTest::TestSuite
{
 private:
  Teuchos::RCP<LineDefinition> lineDefinition;

 public:
  void setUp() override { lineDefinition = Teuchos::rcp(new LineDefinition()); }

  void tearDown() override { lineDefinition = Teuchos::null; }

  // Tag
  void Test_AddTag()
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_TagRequired_But_NothingGiven()
  {
    std::istringstream input("");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_TagRequired_But_IntGiven()
  {
    std::istringstream input("1");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_TagRequired_But_DoubleGiven()
  {
    std::istringstream input("1.23");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // String
  void Test_AddString()
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddString("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_StringRequired_But_NothingGiven()
  {
    std::istringstream input("");
    lineDefinition->AddString("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // Int
  void Test_AddInt()
  {
    std::istringstream input("1");
    lineDefinition->AddInt("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_IntRequired_But_NothingGiven()
  {
    std::istringstream input("");
    lineDefinition->AddInt("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_IntRequired_But_DoubleGiven()
  {
    std::istringstream input("1.23");
    lineDefinition->AddInt("OMEGA");
    AssertThrowsError(input);
  }

  void Test_ReadFalse_When_IntConcatenatedWithSomething()
  {
    std::istringstream input("1*8e+2");
    lineDefinition->AddInt("OMEGA");
    AssertThrowsError(input);
  }

  // Int Vector
  void Test_AddIntVector()
  {
    std::istringstream input("1 2 3");
    lineDefinition->AddIntVector("OMEGA", 3);
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_IntVectorRequired_But_TooFewVectorEntriesGiven()
  {
    std::istringstream input("1 2");
    lineDefinition->AddIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_IntVectorRequired_But_TooManyVectorEntriesGiven()
  {
    std::istringstream input("1 2 3 4");
    lineDefinition->AddIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_IntVectorRequired_But_DoubleVectorEntriesGiven()
  {
    std::istringstream input("1.23 2.34 3.45");
    lineDefinition->AddIntVector("OMEGA", 3);
    AssertThrowsError(input);
  }

  // Double Vector
  void Test_AddDoubleVector()
  {
    std::istringstream input("1 2 3");
    lineDefinition->AddDoubleVector("OMEGA", 3);
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_DoubleVectorRequired_But_TooFewVectorEntriesGiven()
  {
    std::istringstream input("1 2");
    lineDefinition->AddDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_DoubleVectorRequired_But_TooManyVectorEntriesGiven()
  {
    std::istringstream input("1 2 3 4");
    lineDefinition->AddDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  // Named String
  void Test_AddNamedString()
  {
    std::istringstream input("OMEGA TEST");
    lineDefinition->AddNamedString("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_NamedStringRequired_But_NothingGiven()
  {
    std::istringstream input("OMEGA ");
    lineDefinition->AddNamedString("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_NamedStringRequired_But_NoNameGiven()
  {
    std::istringstream input("TEST");
    lineDefinition->AddNamedString("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // Named Int
  void Test_AddNamedInt()
  {
    std::istringstream input("OMEGA 1");
    lineDefinition->AddNamedInt("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_NamedIntRequired_But_NothingGiven()
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddNamedInt("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_NamedIntRequired_But_DoubleGiven()
  {
    std::istringstream input("OMEGA 1.23");
    lineDefinition->AddNamedInt("OMEGA");
    AssertThrowsError(input);
  }

  void Test_ReadFalse_When_NamedIntRequired_But_NoNameGiven()
  {
    std::istringstream input("1");
    lineDefinition->AddNamedInt("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // Named Int Vector
  void Test_AddNamedIntVector()
  {
    std::istringstream input("OMEGA 1 2 3");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_NamedIntVectorRequired_But_TooFewVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1 2");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_NamedIntVectorRequired_But_TooManyVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1 2 3 4");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_NamedIntVectorRequired_But_DoubleVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertThrowsError(input);
  }

  void Test_ReadFalse_When_NamedIntVectorRequired_But_NoNameGiven()
  {
    std::istringstream input("1.23 2.34 3.45");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  // Named Double
  void Test_AddNamedDouble()
  {
    std::istringstream input("OMEGA 1.23");
    lineDefinition->AddNamedDouble("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_NamedDoubleRequired_But_NothingGiven()
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddNamedDouble("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_NamedDoubleConcatenatedWithSomething()
  {
    std::istringstream input("OMEGA 123.45*893");
    lineDefinition->AddNamedDouble("OMEGA");
    AssertThrowsError(input);
  }

  void Test_ReadFalse_When_NamedDoubleRequired_But_NoNameGiven()
  {
    std::istringstream input("1.23");
    lineDefinition->AddNamedDouble("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // Named Double Vector
  void Test_AddNamedDoubleVector()
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    AssertReadReturnsTrue(input);
  }

  void Test_ReadFalse_When_NamedDoubleVectorRequired_But_TooFewVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1.23 2.34");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_NamedDoubleVectorRequired_But_TooManyVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45 4.56");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void Test_ReadFalse_When_NamedDoubleVectorRequired_But_NoNameGiven()
  {
    std::istringstream input("1.23 2.34 3.45");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  // Too much in stream except from comments and whitespaces
  void Test_ReadFalse_When_StreamHasTooMuch()
  {
    std::istringstream input("OMEGA More data that is not in lineDefinition.");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void Test_ReadTrue_When_StreamHasTooMuchButOnlyComments()
  {
    std::istringstream input("OMEGA // Comment");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void Test_ReadTrue_When_StreamHasTooMuchButOnlyWhitespaces()
  {
    std::istringstream input("OMEGA           ");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsTrue(input);
  }

  // Empty stream
  void Test_ReadTrue_When_EmptyStreamIntoEmptyDefinition()
  {
    std::istringstream input("");
    AssertReadReturnsTrue(input);
  }

 private:
  void AssertReadReturnsTrue(std::istringstream& input) { TS_ASSERT(lineDefinition->Read(input)); }

  void AssertReadReturnsFalse(std::istringstream& input)
  {
    TS_ASSERT(not lineDefinition->Read(input));
  }

  void AssertThrowsError(std::istringstream& input)
  {
    TS_ASSERT_THROWS_ANYTHING(lineDefinition->Read(input));
  }
};


#endif