/*! \file

\brief unit tests for LineDefinition

\level 3

\maintainer Sebastian Proell

*/
#ifndef UNIT_DRT_LINEDEFINITION_H_
#define UNIT_DRT_LINEDEFINITION_H_

#include <cxxtest/TestSuite.h>
#include "../../../src/drt_lib/drt_linedefinition.H"

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

  void test_AddTag()
  {
    lineDefinition->AddTag("OMEGA");
    std::istringstream input("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void test_AddString()
  {
    lineDefinition->AddString("OMEGA");
    std::istringstream input("OMEGA");
    AssertReadReturnsTrue(input);
  }

  void test_AddNamedString()
  {
    lineDefinition->AddNamedString("OMEGA");
    std::istringstream input("OMEGA TEST");
    AssertReadReturnsTrue(input);
  }

  void test_AddNamedInt()
  {
    lineDefinition->AddNamedInt("OMEGA");
    std::istringstream input("OMEGA 1");
    AssertReadReturnsTrue(input);
  }

  void test_ReadFalse_When_IntRequiredButDoubleGiven()
  {
    lineDefinition->AddNamedInt("OMEGA");
    std::istringstream input("OMEGA 1.23");
    AssertReadReturnsFalse(input);
  }

  void test_ReadFalse_When_IntRequiredButNothingGiven()
  {
    lineDefinition->AddNamedInt("OMEGA");
    std::istringstream input("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void test_ReadTrue_When_EmptyStreamIntoEmptyDefinition()
  {
    std::istringstream input("");
    AssertReadReturnsTrue(input);
  }

  void test_ReadFalse_When_EmptyStreamIntoFilledDefinition()
  {
    std::istringstream input("");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void test_ReadFalse_When_DoubleConcatenatedWithSomething()
  {
    std::istringstream input("OMEGA 123.45*893");
    lineDefinition->AddNamedDouble("OMEGA");
    AssertReadReturnsFalse(input);
  }

  void test_ReadFalse_When_StreamHasTooMuch()
  {
    std::istringstream input("OMEGA More data that is not in lineDefinition.");
    lineDefinition->AddTag("OMEGA");
    AssertReadReturnsFalse(input);
  }

 private:
  void AssertReadReturnsTrue(std::istringstream& input) { TS_ASSERT(lineDefinition->Read(input)); }
  void AssertReadReturnsFalse(std::istringstream& input)
  {
    TS_ASSERT(not lineDefinition->Read(input));
  }
};


#endif
