/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for FunctionVariable classes
\level 3
*/


#ifndef UNIT_DRT_FUNCTION_VARIABLES_H
#define UNIT_DRT_FUNCTION_VARIABLES_H

#include "Unittests/src/common/unit_cxx_test_wrapper.H"
#include "../../../src/drt_lib/drt_functionvariables.H"

namespace DRT::UTILS
{
  class PiecewiseVariable_TestSuite;
}  // namespace DRT::UTILS



class DRT::UTILS::PiecewiseVariable_TestSuite : public BACICxxTestWrapper
{
 private:
  class FunctionVariableStub : public FunctionVariable
  {
   public:
    FunctionVariableStub(double start, double end, double returnValue = 0.0)
        : FunctionVariable("stub"), start(start), end(end), returnValue(returnValue)
    {
    }

    double Value(const double t) override { return returnValue; }

    double TimeDerivativeValue(const double t, const unsigned int deg) override
    {
      return returnValue;
    }

    bool ContainTime(const double t) override { return t >= start && t <= end; }

   private:
    const double start;
    const double end;
    const double returnValue;
  };

 public:
  void TestContainTimeSingle()
  {
    auto stub = Teuchos::rcp(new FunctionVariableStub(0.0, 1.0));
    PiecewiseVariable piecewiseVariable("test", {stub});

    TS_ASSERT(piecewiseVariable.ContainTime(0.0));
    TS_ASSERT(piecewiseVariable.ContainTime(0.5));
    TS_ASSERT(piecewiseVariable.ContainTime(1.0));
    TS_ASSERT(not piecewiseVariable.ContainTime(2.0));
  }

  void TestValueSingle()
  {
    auto stub = Teuchos::rcp(new FunctionVariableStub(0.0, 1.0));
    PiecewiseVariable piecewiseVariable("test", {stub});

    TS_ASSERT_EQUALS(piecewiseVariable.Value(0.0), 0.0);
    TS_ASSERT_EQUALS(piecewiseVariable.Value(0.5), 0.0);
    TS_ASSERT_EQUALS(piecewiseVariable.Value(1.0), 0.0);
    TS_ASSERT_THROWS_ANYTHING(piecewiseVariable.Value(2.0));
  }

  void TestTimeDerivativeValueSingle()
  {
    auto stub = Teuchos::rcp(new FunctionVariableStub(0.0, 1.0));
    PiecewiseVariable piecewiseVariable("test", {stub});

    TS_ASSERT_EQUALS(piecewiseVariable.TimeDerivativeValue(0.0, 1), 0.0);
    TS_ASSERT_EQUALS(piecewiseVariable.TimeDerivativeValue(0.5, 1), 0.0);
    TS_ASSERT_EQUALS(piecewiseVariable.TimeDerivativeValue(1.0, 1), 0.0);
    TS_ASSERT_THROWS_ANYTHING(piecewiseVariable.TimeDerivativeValue(2.0, 1));
  }

  void TestContainTimeMutliple()
  {
    // define a piece-wise variable which has overlapping areas
    PiecewiseVariable piecewiseVariable(
        "test", {Teuchos::rcp(new FunctionVariableStub(0.0, 1.0, 1.0)),
                    Teuchos::rcp(new FunctionVariableStub(0.5, 2.0, 2.0)),
                    Teuchos::rcp(new FunctionVariableStub(1.0, 3.0, 3.0))});

    TS_ASSERT(piecewiseVariable.ContainTime(0.0));
    TS_ASSERT(piecewiseVariable.ContainTime(1.0));
    TS_ASSERT(piecewiseVariable.ContainTime(2.0));
    TS_ASSERT(piecewiseVariable.ContainTime(3.0));

    TS_ASSERT(not piecewiseVariable.ContainTime(4.0));
  }


  void TestValueMultiple()
  {
    // define a piece-wise variable which has overlapping areas
    PiecewiseVariable piecewiseVariable(
        "test", {Teuchos::rcp(new FunctionVariableStub(0.0, 1.0, 1.0)),
                    Teuchos::rcp(new FunctionVariableStub(0.5, 2.0, 2.0)),
                    Teuchos::rcp(new FunctionVariableStub(1.0, 3.0, 3.0))});

    TS_ASSERT_EQUALS(piecewiseVariable.Value(0.0), 1.0);
    TS_ASSERT_EQUALS(piecewiseVariable.Value(0.8), 1.0);
    TS_ASSERT_EQUALS(piecewiseVariable.Value(1.5), 2.0);
    TS_ASSERT_EQUALS(piecewiseVariable.Value(2.1), 3.0);
  }

  void TestTimeDerivativeValueMultiple()
  {
    // define a piece-wise variable which has overlapping areas
    PiecewiseVariable piecewiseVariable(
        "test", {Teuchos::rcp(new FunctionVariableStub(0.0, 1.0, 1.0)),
                    Teuchos::rcp(new FunctionVariableStub(0.5, 2.0, 2.0)),
                    Teuchos::rcp(new FunctionVariableStub(1.0, 3.0, 3.0))});

    TS_ASSERT_EQUALS(piecewiseVariable.TimeDerivativeValue(0.0, 0), 1.0);
    TS_ASSERT_EQUALS(piecewiseVariable.TimeDerivativeValue(0.8, 0), 1.0);
    TS_ASSERT_EQUALS(piecewiseVariable.TimeDerivativeValue(1.5, 0), 2.0);
    TS_ASSERT_EQUALS(piecewiseVariable.TimeDerivativeValue(2.1, 0), 3.0);
  }
};


#endif
