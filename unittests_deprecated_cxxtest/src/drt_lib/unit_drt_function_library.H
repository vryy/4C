/*---------------------------------------------------------------------------*/
/*! \file

\brief unittests for the function library

\level 3

*/


#ifndef UNIT_DRT_FUNCTION_LIBRARY_H
#define UNIT_DRT_FUNCTION_LIBRARY_H

#include "src/common/unit_cxx_test_wrapper.H"
#include "../../../src/drt_lib/drt_function.H"
#include "../../../src/drt_lib/drt_function_library.H"

namespace DRT
{
  namespace UTILS
  {
    class TranslatedFunction_TestSuite;
  }
}  // namespace DRT



class DRT::UTILS::TranslatedFunction_TestSuite : public BACICxxTestWrapper
{
 private:
  constexpr static double A_CONSTANT = 1.0;

  Teuchos::RCP<Function> constant_;
  Teuchos::RCP<Function> scalarField_;
  Teuchos::RCP<Function> linearTime_;

  const std::array<double, 3> at000_ = {0.0, 0.0, 0.0};
  const std::array<double, 3> at111_ = {1.0, 1.0, 1.0};
  const std::array<double, 3> at123_ = {1.0, 2.0, 3.0};

  /// a time- and space-independent constant vector function
  class ConstantVectorFunctionStub : public DRT::UTILS::Function
  {
   public:
    double Evaluate(const int index, const double* x, double t) override { return A_CONSTANT; }
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, double t, const unsigned int deg) override
    {
      return {0.0, 0.0};
    }
    std::vector<double> EvaluateSpatialDerivative(
        const int index, const double* x, const double t) override
    {
      return {0.0, 0.0, 0.0};
    }
    std::size_t NumberComponents() override { return 3; }
  };

  /// a scalar field f(x) = x_1^2+2*x_2^2+4*x_3^2
  class ScalarFieldFunctionStub : public DRT::UTILS::Function
  {
   public:
    double Evaluate(const int index, const double* x, double t) override
    {
      return x[0] * x[0] + 2 * x[1] * x[1] + 4 * x[2] * x[2];
    }
    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, double t, const unsigned int deg) override
    {
      auto result = std::vector<double>(2, 0.0);
      result[0] = Evaluate(index, x, t);
      return result;
    }
    std::vector<double> EvaluateSpatialDerivative(
        const int index, const double* x, const double t) override
    {
      return {2 * x[0], 4 * x[1], 8 * x[2]};
    }
    std::size_t NumberComponents() override { return 1; }
  };

  /// a time dependent vector function f(t) = [t, t, t]
  class LinearTimeVectorFunctionStub : public DRT::UTILS::Function
  {
   public:
    double Evaluate(const int index, const double* x, double t) override { return t; }

    std::vector<double> EvaluateTimeDerivative(
        const int index, const double* x, double t, const unsigned int deg) override
    {
      auto result = std::vector<double>(2, 0.0);
      result[0] = Evaluate(index, x, t);
      result[1] = 1;
      return result;
    }
    std::vector<double> EvaluateSpatialDerivative(
        const int index, const double* x, const double t) override
    {
      return {0.0, 0.0, 0.0};
    }
    std::size_t NumberComponents() override { return 3; }
  };


 public:
  void Setup() override
  {
    constant_ = Teuchos::rcp(new ConstantVectorFunctionStub());
    scalarField_ = Teuchos::rcp(new ScalarFieldFunctionStub());
    linearTime_ = Teuchos::rcp(new LinearTimeVectorFunctionStub());
  }

  void TearDown() override
  {
    constant_ = Teuchos::null;
    linearTime_ = Teuchos::null;
    scalarField_ = Teuchos::null;
  }

  void TestEvaluateConstantOriginConstantLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(constant_, constant_));
    // the whole function is a constant in all components
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111_.data(), 0), A_CONSTANT, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(1, at123_.data(), 0), A_CONSTANT, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(2, at111_.data(), 0), A_CONSTANT, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123_.data(), 1000), A_CONSTANT, 1.0e-15);
  }

  void TestEvaluateConstantOriginSpaceDependentLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(constant_, scalarField_));
    // the local origin is fixed at 1,1,1 and does not move
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111_.data(), 0), 0, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111_.data(), 1), 0, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123_.data(), 0),
        Expected(at123_.data(), WithOrigin(at111_.data())), 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123_.data(), 1),
        Expected(at123_.data(), WithOrigin(at111_.data())), 1.0e-15);
  }

  void TestEvaluateMovingOriginConstantLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime_, constant_));
    // local function is a constant and origin does not matter
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111_.data(), 0), A_CONSTANT, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111_.data(), 1), A_CONSTANT, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123_.data(), 0), A_CONSTANT, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123_.data(), 1), A_CONSTANT, 1.0e-15);
  }

  void TestEvaluateMovingOriginSpaceDependentLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime_, scalarField_));
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at111_.data(), 0),
        Expected(at111_.data(), WithOrigin(at000_.data())), 1.0e-15);
    TS_ASSERT_DELTA(testFunction->Evaluate(0, at123_.data(), 1),
        Expected(at123_.data(), WithOrigin(at111_.data())), 1.0e-15);
  }

  void TestEvaluateTimeDerivativeConstantOriginConstantLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(constant_, constant_));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 1)[1], 0, 1.0e-15);
  }

  void TestEvaluateTimeDerivativeConstantOriginSpaceDependentLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(constant_, scalarField_));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 1)[1], 0, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123_.data(), 0, 1)[1], 0, 1.0e-15);
  }

  void TestEvaluateTimeDerivativeMovingOriginConstantLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime_, constant_));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 1)[1], 0, 1.0e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123_.data(), 1, 1)[1], 0, 1.0e-15);
  }

  void TestEvaluateTimeDerivativeMovingOriginSpaceDependentLocal()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime_, scalarField_));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 1)[1],
        ExpectedTimeDerivative(at111_.data(), 0), 1.0e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at111_.data(), 1, 1)[1],
        ExpectedTimeDerivative(at111_.data(), 1), 1.0e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123_.data(), 0, 1)[1],
        ExpectedTimeDerivative(at123_.data(), 0), 1.0e-15);
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123_.data(), 1, 1)[1],
        ExpectedTimeDerivative(at123_.data(), 1), 1.0e-15);
  }

  void TestEvaluateTimeDerivativeDegree0()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime_, scalarField_));
    TS_ASSERT_DELTA(testFunction->EvaluateTimeDerivative(0, at123_.data(), 1, 0)[0],
        Expected(at123_.data(), WithOrigin(at111_.data())), 1.0e-15);
  }

  void TestShouldThrowWhenEvaluateTimeDerivativeDegreeGreater1()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime_, scalarField_));
    TS_ASSERT_THROWS_ANYTHING(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 2));
  }

  void TestShouldThrowWhenEvaluatedComponentNegative()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime_, scalarField_));
    TS_ASSERT_THROWS_ANYTHING(testFunction->Evaluate(-1, at111_.data(), 0));
  }

  void TestShouldThrowWhenEvaluatedComponentTooLarge()
  {
    Teuchos::RCP<Function> testFunction =
        Teuchos::rcp(new TranslatedFunction(linearTime_, scalarField_));
    TS_ASSERT_THROWS_ANYTHING(testFunction->Evaluate(2, at111_.data(), 0));
  }

  void TestShouldThrowWhenOriginHasOnlyOneComponent()
  {
    TS_ASSERT_THROWS_ANYTHING(new TranslatedFunction(scalarField_, constant_));
  }

 private:
  inline double Expected(const double* at, const double* withOrigin)
  {
    return (at[0] - withOrigin[0]) * (at[0] - withOrigin[0]) +
           2 * (at[1] - withOrigin[1]) * (at[1] - withOrigin[1]) +
           4 * (at[2] - withOrigin[2]) * (at[2] - withOrigin[2]);
  }

  inline double ExpectedTimeDerivative(const double* at, const double time)
  {
    return -2 * (at[0] - time) - 4 * (at[1] - time) - 8 * (at[2] - time);
  }

  /// for readability
  inline const double* WithOrigin(const double* at) { return at; }
};


#endif
