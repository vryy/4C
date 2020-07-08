/*---------------------------------------------------------------------------*/
/*! \file

\brief unittests for the function library

\level 3

*/


#ifndef UNIT_DRT_FUNCTION_LIBRARY_H
#define UNIT_DRT_FUNCTION_LIBRARY_H

#include <cxxtest/TestSuite.h>
#include "../../../src/drt_lib/drt_function.H"
#include "../../../src/drt_lib/drt_function_library.H"

namespace DRT
{
  namespace UTILS
  {
    class TranslatedFunction_TestSuite;
  }
}  // namespace DRT



class DRT::UTILS::TranslatedFunction_TestSuite : public CxxTest::TestSuite
{
 private:
  constexpr static double A_CONSTANT = 1.0;

  Teuchos::RCP<Function> constant;
  Teuchos::RCP<Function> scalarField;
  Teuchos::RCP<Function> linearTime;

  const double at000[3] = {0, 0, 0};
  const double at111[3] = {1, 1, 1};
  const double at123[3] = {1, 2, 3};

  /// a time- and space-independent constant vector function
  class ConstantVectorFunctionStub : public DRT::UTILS::Function
  {
   public:
    double Evaluate(const int index, const double* x, double t) override { return A_CONSTANT; };
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, double t, const unsigned int deg) override
    {
      return std::vector<double>(2, 0.0);
    };
    std::vector<double> EvaluateSpatialDerivative(
        const int index, const double* x, const double t) override
    {
      return std::vector<double>(3, 0);
    }
    int NumberComponents() override { return 3; };
  };

  /// a scalar field f(x) = x_1^2+2*x_2^2+4*x_3^2
  class ScalarFieldFunctionStub : public DRT::UTILS::Function
  {
   public:
    double Evaluate(const int index, const double* x, double t) override
    {
      return (x[0] * x[0] + 2 * x[1] * x[1] + 4 * x[2] * x[2]);
    };
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, double t, const unsigned int deg) override
    {
      auto result = std::vector<double>(2, 0.0);
      result[0] = Evaluate(index, x, t);
      return result;
    };
    std::vector<double> EvaluateSpatialDerivative(
        const int index, const double* x, const double t) override
    {
      return std::vector<double>{2 * x[0], 4 * x[1], 8 * x[2]};
    }
    int NumberComponents() override { return 1; };
  };

  /// a time dependent vector function f(t) = [t, t, t]
  class LinearTimeVectorFunctionStub : public DRT::UTILS::Function
  {
   public:
    double Evaluate(const int index, const double* x, double t) override { return t; };

    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, double t, const unsigned int deg) override
    {
      auto result = std::vector<double>(2, 0.0);
      result[0] = Evaluate(index, x, t);
      result[1] = 1;
      return result;
    };
    std::vector<double> EvaluateSpatialDerivative(
        const int index, const double* x, const double t) override
    {
      return std::vector<double>(3, 0);
    }
    int NumberComponents() override { return 3; };
  };


 public:
  void setUp() override
  {
    constant = Teuchos::rcp(new ConstantVectorFunctionStub());
    scalarField = Teuchos::rcp(new ScalarFieldFunctionStub());
    linearTime = Teuchos::rcp(new LinearTimeVectorFunctionStub());
  };

  void tearDown() override
  {
    constant = Teuchos::null;
    linearTime = Teuchos::null;
    scalarField = Teuchos::null;
  };

  void test_Evaluate_ConstantOrigin_ConstantLocal()
  {
    Teuchos::RCP<Function> testFunction = Teuchos::rcp(new TranslatedFunction(constant, constant));
    // the whole function is a constant in all components
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111, 0), A_CONSTANT, 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(1, at123, 0), A_CONSTANT, 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(2, at111, 0), A_CONSTANT, 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123, 1000), A_CONSTANT, 1e-15);
  }

  void test_Evaluate_ConstantOrigin_SpaceDependentLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(constant, scalarField));
    // the local origin is fixed at 1,1,1 and does not move
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111, 0), 0, 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111, 1), 0, 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123, 0), expected(at123, withOrigin(at111)), 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123, 1), expected(at123, withOrigin(at111)), 1e-15);
  }

  void test_Evaluate_MovingOrigin_ConstantLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime, constant));
    // local function is a constant and origin does not matter
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111, 0), A_CONSTANT, 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111, 1), A_CONSTANT, 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123, 0), A_CONSTANT, 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123, 1), A_CONSTANT, 1e-15);
  }

  void test_Evaluate_MovingOrigin_SpaceDependentLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime, scalarField));
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111, 0), expected(at111, withOrigin(at000)), 1e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123, 1), expected(at123, withOrigin(at111)), 1e-15);
  }

  void test_EvaluateTimeDerivative_ConstantOrigin_ConstantLocal()
  {
    Teuchos::RCP<Function> testFunction = Teuchos::rcp(new TranslatedFunction(constant, constant));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111, 0, 1)[1], 0, 1e-15);
  }

  void test_EvaluateTimeDerivative_ConstantOrigin_SpaceDependentLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(constant, scalarField));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111, 0, 1)[1], 0, 1e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123, 0, 1)[1], 0, 1e-15);
  }

  void test_EvaluateTimeDerivative_MovingOrigin_ConstantLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime, constant));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111, 0, 1)[1], 0, 1e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123, 1, 1)[1], 0, 1e-15);
  }

  void test_EvaluateTimeDerivative_MovingOrigin_SpaceDependentLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime, scalarField));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111, 0, 1)[1],
        expectedTimeDerivative(at111, 0), 1e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111, 1, 1)[1],
        expectedTimeDerivative(at111, 1), 1e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123, 0, 1)[1],
        expectedTimeDerivative(at123, 0), 1e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123, 1, 1)[1],
        expectedTimeDerivative(at123, 1), 1e-15);
  }

  void test_EvaluateTimeDerivative_Degree0()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime, scalarField));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123, 1, 0)[0],
        expected(at123, withOrigin(at111)), 1e-15);
  }

  void test_ShouldThrow_When_EvaluateTimeDerivative_DegreeGreater1()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime, scalarField));
    TS_ASSERT_THROWS_ANYTHING(testFunction->EvaluateTimeDerivative(0, at111, 0, 2));
  }

  void test_ShouldThrow_When_EvaluatedComponentNegative()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime, scalarField));
    TS_ASSERT_THROWS_ANYTHING(testFunction->Evaluate(-1, at111, 0));
  }

  void test_ShouldThrow_When_EvaluatedComponentTooLarge()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime, scalarField));
    TS_ASSERT_THROWS_ANYTHING(testFunction->Evaluate(2, at111, 0));
  }

  void test_ShouldThrow_When_OriginHasOnlyOneComponent()
  {
    TS_ASSERT_THROWS_ANYTHING(new TranslatedFunction(scalarField, constant));
  }

 private:
  inline double expected(const double* at, const double* withOrigin)
  {
    return (at[0] - withOrigin[0]) * (at[0] - withOrigin[0]) +
           2 * (at[1] - withOrigin[1]) * (at[1] - withOrigin[1]) +
           4 * (at[2] - withOrigin[2]) * (at[2] - withOrigin[2]);
  }

  inline double expectedTimeDerivative(const double* at, const double time)
  {
    return -2 * (at[0] - time) - 4 * (at[1] - time) - 8 * (at[2] - time);
  }

  /// for readability
  inline const double* withOrigin(const double* at) { return at; }
};


#endif
