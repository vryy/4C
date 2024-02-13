/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for FunctionVariable classes

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_utils_functionvariables.hpp"

BACI_NAMESPACE_OPEN

namespace
{
  class FunctionVariableStub : public CORE::UTILS::FunctionVariable
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

  TEST(PiecewiseVariableTest, ContainTimeSingle)
  {
    auto stub = Teuchos::rcp(new FunctionVariableStub(0.0, 1.0));
    CORE::UTILS::PiecewiseVariable piecewiseVariable("test", {stub});

    EXPECT_TRUE(piecewiseVariable.ContainTime(0.0));
    EXPECT_TRUE(piecewiseVariable.ContainTime(0.5));
    EXPECT_TRUE(piecewiseVariable.ContainTime(1.0));
    EXPECT_TRUE(not piecewiseVariable.ContainTime(2.0));
  }

  TEST(PiecewiseVariableTest, ValueSingle)
  {
    auto stub = Teuchos::rcp(new FunctionVariableStub(0.0, 1.0));
    CORE::UTILS::PiecewiseVariable piecewiseVariable("test", {stub});

    EXPECT_EQ(piecewiseVariable.Value(0.0), 0.0);
    EXPECT_EQ(piecewiseVariable.Value(0.5), 0.0);
    EXPECT_EQ(piecewiseVariable.Value(1.0), 0.0);
    EXPECT_ANY_THROW(piecewiseVariable.Value(2.0));
  }

  TEST(PiecewiseVariableTest, TimeDerivativeValueSingle)
  {
    auto stub = Teuchos::rcp(new FunctionVariableStub(0.0, 1.0));
    CORE::UTILS::PiecewiseVariable piecewiseVariable("test", {stub});

    EXPECT_EQ(piecewiseVariable.TimeDerivativeValue(0.0, 1), 0.0);
    EXPECT_EQ(piecewiseVariable.TimeDerivativeValue(0.5, 1), 0.0);
    EXPECT_EQ(piecewiseVariable.TimeDerivativeValue(1.0, 1), 0.0);
    EXPECT_ANY_THROW(piecewiseVariable.TimeDerivativeValue(2.0, 1));
  }

  TEST(PiecewiseVariableTest, ContainTimeMultiple)
  {
    // define a piece-wise variable which has overlapping areas
    CORE::UTILS::PiecewiseVariable piecewiseVariable(
        "test", {Teuchos::rcp(new FunctionVariableStub(0.0, 1.0, 1.0)),
                    Teuchos::rcp(new FunctionVariableStub(0.5, 2.0, 2.0)),
                    Teuchos::rcp(new FunctionVariableStub(1.0, 3.0, 3.0))});

    EXPECT_TRUE(piecewiseVariable.ContainTime(0.0));
    EXPECT_TRUE(piecewiseVariable.ContainTime(1.0));
    EXPECT_TRUE(piecewiseVariable.ContainTime(2.0));
    EXPECT_TRUE(piecewiseVariable.ContainTime(3.0));

    EXPECT_TRUE(not piecewiseVariable.ContainTime(4.0));
  }

  TEST(PiecewiseVariableTest, ValueMultiple)
  {
    // define a piece-wise variable which has overlapping areas
    CORE::UTILS::PiecewiseVariable piecewiseVariable(
        "test", {Teuchos::rcp(new FunctionVariableStub(0.0, 1.0, 1.0)),
                    Teuchos::rcp(new FunctionVariableStub(0.5, 2.0, 2.0)),
                    Teuchos::rcp(new FunctionVariableStub(1.0, 3.0, 3.0))});

    EXPECT_EQ(piecewiseVariable.Value(0.0), 1.0);
    EXPECT_EQ(piecewiseVariable.Value(0.8), 1.0);
    EXPECT_EQ(piecewiseVariable.Value(1.5), 2.0);
    EXPECT_EQ(piecewiseVariable.Value(2.1), 3.0);
  }

  TEST(PiecewiseVariableTest, TimeDerivativeValueMultiple)
  {
    // define a piece-wise variable which has overlapping areas
    CORE::UTILS::PiecewiseVariable piecewiseVariable(
        "test", {Teuchos::rcp(new FunctionVariableStub(0.0, 1.0, 1.0)),
                    Teuchos::rcp(new FunctionVariableStub(0.5, 2.0, 2.0)),
                    Teuchos::rcp(new FunctionVariableStub(1.0, 3.0, 3.0))});

    EXPECT_EQ(piecewiseVariable.TimeDerivativeValue(0.0, 0), 1.0);
    EXPECT_EQ(piecewiseVariable.TimeDerivativeValue(0.8, 0), 1.0);
    EXPECT_EQ(piecewiseVariable.TimeDerivativeValue(1.5, 0), 2.0);
    EXPECT_EQ(piecewiseVariable.TimeDerivativeValue(2.1, 0), 3.0);
  }
}  // namespace

BACI_NAMESPACE_CLOSE
