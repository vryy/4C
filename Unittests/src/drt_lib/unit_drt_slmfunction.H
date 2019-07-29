/*---------------------------------------------------------------------------*/
/*!

\brief unittests for slm functions

\level 3

\maintainer  Sebastian Proell
*/


#ifndef UNIT_DRT_SLMFUNCTION_H
#define UNIT_DRT_SLMFUNCTION_H

#include <cxxtest/TestSuite.h>
#include "../../../src/drt_lib/drt_function.H"
#include "../../../src/drt_lib/drt_slmfunction.H"

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
    int NumberComponents() override { return 1; };
  };

  /// a time dependent vector function f(t) = [t, t, t]
  class LinearTimeVectorFunctionStub : public DRT::UTILS::Function
  {
   public:
    double Evaluate(const int index, const double* x, double t) override { return t; };
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
  double expected(const double* at, const double* withOrigin)
  {
    return (at[0] - withOrigin[0]) * (at[0] - withOrigin[0]) +
           2 * (at[1] - withOrigin[1]) * (at[1] - withOrigin[1]) +
           4 * (at[2] - withOrigin[2]) * (at[2] - withOrigin[2]);
  }

  /// for readability
  const double* withOrigin(const double* at) { return at; }
};


#endif
