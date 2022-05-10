/*! \file

\brief unit tests for LineDefinition

\level 3
*/
#ifndef UNIT_DRT_LINEDEFINITION_H_
#define UNIT_DRT_LINEDEFINITION_H_

#include "src/common/unit_cxx_test_wrapper.H"
#include "src/drt_lib/drt_linedefinition.H"

namespace DRT
{
  namespace INPUT
  {
    class LineDefinition_TestSuite;
  }
}  // namespace DRT


class DRT::INPUT::LineDefinition_TestSuite : public BACICxxTestWrapper
{
 private:
  Teuchos::RCP<LineDefinition> lineDefinition;

 public:
  void Setup() override { lineDefinition = Teuchos::rcp(new LineDefinition()); }

  void TearDown() override { lineDefinition = Teuchos::null; }

  // Tag
  void TestAddTag()
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenTagRequiredButNothingGiven()
  {
    std::istringstream input("");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenTagRequiredButIntGiven()
  {
    std::istringstream input("1");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenTagRequiredButDoubleGiven()
  {
    std::istringstream input("1.23");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // String
  void TestAddString()
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddString("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenStringRequiredButNothingGiven()
  {
    std::istringstream input("");
    lineDefinition->AddString("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // Int
  void TestAddInt()
  {
    std::istringstream input("1");
    lineDefinition->AddInt("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenIntRequiredButNothingGiven()
  {
    std::istringstream input("");
    lineDefinition->AddInt("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenIntRequiredButDoubleGiven()
  {
    std::istringstream input("1.23");
    lineDefinition->AddInt("OMEGA");
    AssertThrowsError(input);
  }

  void TestReadFalseWhenIntConcatenatedWithSomething()
  {
    std::istringstream input("1*8e+2");
    lineDefinition->AddInt("OMEGA");
    AssertThrowsError(input);
  }

  // Int Vector
  void TestAddIntVector()
  {
    std::istringstream input("1 2 3");
    lineDefinition->AddIntVector("OMEGA", 3);
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenIntVectorRequiredButTooFewVectorEntriesGiven()
  {
    std::istringstream input("1 2");
    lineDefinition->AddIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenIntVectorRequiredButTooManyVectorEntriesGiven()
  {
    std::istringstream input("1 2 3 4");
    lineDefinition->AddIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenIntVectorRequiredButDoubleVectorEntriesGiven()
  {
    std::istringstream input("1.23 2.34 3.45");
    lineDefinition->AddIntVector("OMEGA", 3);
    AssertThrowsError(input);
  }

  // Double Vector
  void TestAddDoubleVector()
  {
    std::istringstream input("1 2 3");
    lineDefinition->AddDoubleVector("OMEGA", 3);
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenDoubleVectorRequiredButTooFewVectorEntriesGiven()
  {
    std::istringstream input("1 2");
    lineDefinition->AddDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenDoubleVectorRequiredButTooManyVectorEntriesGiven()
  {
    std::istringstream input("1 2 3 4");
    lineDefinition->AddDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  // Named String
  void TestAddNamedString()
  {
    std::istringstream input("OMEGA TEST");
    lineDefinition->AddNamedString("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenNamedStringRequiredButNothingGiven()
  {
    std::istringstream input("OMEGA ");
    lineDefinition->AddNamedString("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenNamedStringRequiredButNoNameGiven()
  {
    std::istringstream input("TEST");
    lineDefinition->AddNamedString("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // Named Int
  void TestAddNamedInt()
  {
    std::istringstream input("OMEGA 1");
    lineDefinition->AddNamedInt("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenNamedIntRequiredButNothingGiven()
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddNamedInt("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenNamedIntRequiredButDoubleGiven()
  {
    std::istringstream input("OMEGA 1.23");
    lineDefinition->AddNamedInt("OMEGA");
    AssertThrowsError(input);
  }

  void TestReadFalseWhenNamedIntRequiredButNoNameGiven()
  {
    std::istringstream input("1");
    lineDefinition->AddNamedInt("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // Named Int Vector
  void TestAddNamedIntVector()
  {
    std::istringstream input("OMEGA 1 2 3");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenNamedIntVectorRequiredButTooFewVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1 2");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenNamedIntVectorRequiredButTooManyVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1 2 3 4");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenNamedIntVectorRequiredButDoubleVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertThrowsError(input);
  }

  void TestReadFalseWhenNamedIntVectorRequiredButNoNameGiven()
  {
    std::istringstream input("1.23 2.34 3.45");
    lineDefinition->AddNamedIntVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  // Named Double
  void TestAddNamedDouble()
  {
    std::istringstream input("OMEGA 1.23");
    lineDefinition->AddNamedDouble("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenNamedDoubleRequiredButNothingGiven()
  {
    std::istringstream input("OMEGA");
    lineDefinition->AddNamedDouble("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenNamedDoubleConcatenatedWithSomething()
  {
    std::istringstream input("OMEGA 123.45*893");
    lineDefinition->AddNamedDouble("OMEGA");
    AssertThrowsError(input);
  }

  void TestReadFalseWhenNamedDoubleRequiredButNoNameGiven()
  {
    std::istringstream input("1.23");
    lineDefinition->AddNamedDouble("OMEGA");
    AssertReadReturnsFalse(input);
  }

  // Named Double Vector
  void TestAddNamedDoubleVector()
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    AssertReadReturnsTrue(input);
  }

  void TestReadFalseWhenNamedDoubleVectorRequiredButTooFewVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1.23 2.34");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenNamedDoubleVectorRequiredButTooManyVectorEntriesGiven()
  {
    std::istringstream input("OMEGA 1.23 2.34 3.45 4.56");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  void TestReadFalseWhenNamedDoubleVectorRequiredButNoNameGiven()
  {
    std::istringstream input("1.23 2.34 3.45");
    lineDefinition->AddNamedDoubleVector("OMEGA", 3);
    AssertReadReturnsFalse(input);
  }

  // Too much in stream except from comments and whitespaces
  void TestReadFalseWhenStreamHasTooMuch()
  {
    std::istringstream input("OMEGA More data that is not in lineDefinition.");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void TestReadTrueWhenStreamHasTooMuchButOnlyComments()
  {
    std::istringstream input("OMEGA // Comment");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void TestReadTrue_When_StreamHasTooMuchButOnlyWhitespaces()
  {
    std::istringstream input("OMEGA           ");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsTrue(input);
  }

  // Empty stream
  void TestReadTrueWhenEmptyStreamIntoEmptyDefinition()
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