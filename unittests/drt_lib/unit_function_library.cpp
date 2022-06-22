/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for the function library

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include "src/drt_lib/drt_function.H"
#include "src/drt_lib/drt_function_library.H"

namespace
{
  class TranslatedFunctionTest : public testing::Test
  {
   public:
    TranslatedFunctionTest()
    {
      constant_ = Teuchos::rcp(new ConstantVectorFunctionStub());
      scalarField_ = Teuchos::rcp(new ScalarFieldFunctionStub());
      linearTime_ = Teuchos::rcp(new LinearTimeVectorFunctionStub());
    }

   protected:
    constexpr static double A_CONSTANT = 1.0;

    Teuchos::RCP<DRT::UTILS::Function> constant_;
    Teuchos::RCP<DRT::UTILS::Function> scalarField_;
    Teuchos::RCP<DRT::UTILS::Function> linearTime_;

    const std::array<double, 3> at000_ = {0.0, 0.0, 0.0};
    const std::array<double, 3> at111_ = {1.0, 1.0, 1.0};
    const std::array<double, 3> at123_ = {1.0, 2.0, 3.0};

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

   private:
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
  };

  TEST_F(TranslatedFunctionTest, EvaluateConstantOriginConstantLocal)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(constant_, constant_));
    // the whole function is a constant in all components
    EXPECT_NEAR(testFunction->Evaluate(0, at111_.data(), 0), A_CONSTANT, 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(1, at123_.data(), 0), A_CONSTANT, 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(2, at111_.data(), 0), A_CONSTANT, 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(0, at123_.data(), 1000), A_CONSTANT, 1.0e-15);
  }

  TEST_F(TranslatedFunctionTest, EvaluateConstantOriginSpaceDependentLocal)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(constant_, scalarField_));
    // the local origin is fixed at 1,1,1 and does not move
    EXPECT_NEAR(testFunction->Evaluate(0, at111_.data(), 0), 0, 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(0, at111_.data(), 1), 0, 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(0, at123_.data(), 0),
        Expected(at123_.data(), WithOrigin(at111_.data())), 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(0, at123_.data(), 1),
        Expected(at123_.data(), WithOrigin(at111_.data())), 1.0e-15);
  }

  TEST_F(TranslatedFunctionTest, EvaluateMovingOriginConstantLocal)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(linearTime_, constant_));
    // local function is a constant and origin does not matter
    EXPECT_NEAR(testFunction->Evaluate(0, at111_.data(), 0), A_CONSTANT, 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(0, at111_.data(), 1), A_CONSTANT, 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(0, at123_.data(), 0), A_CONSTANT, 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(0, at123_.data(), 1), A_CONSTANT, 1.0e-15);
  }

  TEST_F(TranslatedFunctionTest, EvaluateMovingOriginSpaceDependentLocal)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(linearTime_, scalarField_));
    EXPECT_NEAR(testFunction->Evaluate(0, at111_.data(), 0),
        Expected(at111_.data(), WithOrigin(at000_.data())), 1.0e-15);
    EXPECT_NEAR(testFunction->Evaluate(0, at123_.data(), 1),
        Expected(at123_.data(), WithOrigin(at111_.data())), 1.0e-15);
  }

  TEST_F(TranslatedFunctionTest, EvaluateTimeDerivativeConstantOriginConstantLocal)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(constant_, constant_));
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 1)[1], 0, 1.0e-15);
  }

  TEST_F(TranslatedFunctionTest, EvaluateTimeDerivativeConstantOriginSpaceDependentLocal)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(constant_, scalarField_));
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 1)[1], 0, 1.0e-15);
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at123_.data(), 0, 1)[1], 0, 1.0e-15);
  }

  TEST_F(TranslatedFunctionTest, EvaluateTimeDerivativeMovingOriginConstantLocal)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(linearTime_, constant_));
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 1)[1], 0, 1.0e-15);
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at123_.data(), 1, 1)[1], 0, 1.0e-15);
  }

  TEST_F(TranslatedFunctionTest, EvaluateTimeDerivativeMovingOriginSpaceDependentLocal)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(linearTime_, scalarField_));
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 1)[1],
        ExpectedTimeDerivative(at111_.data(), 0), 1.0e-15);
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at111_.data(), 1, 1)[1],
        ExpectedTimeDerivative(at111_.data(), 1), 1.0e-15);
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at123_.data(), 0, 1)[1],
        ExpectedTimeDerivative(at123_.data(), 0), 1.0e-15);
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at123_.data(), 1, 1)[1],
        ExpectedTimeDerivative(at123_.data(), 1), 1.0e-15);
  }

  TEST_F(TranslatedFunctionTest, EvaluateTimeDerivativeDegree0)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(linearTime_, scalarField_));
    EXPECT_NEAR(testFunction->EvaluateTimeDerivative(0, at123_.data(), 1, 0)[0],
        Expected(at123_.data(), WithOrigin(at111_.data())), 1.0e-15);
  }

  TEST_F(TranslatedFunctionTest, ShouldThrowWhenEvaluateTimeDerivativeDegreeGreater1)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(linearTime_, scalarField_));
    EXPECT_ANY_THROW(testFunction->EvaluateTimeDerivative(0, at111_.data(), 0, 2));
  }

  TEST_F(TranslatedFunctionTest, ShouldThrowWhenEvaluatedComponentNegative)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(linearTime_, scalarField_));
    EXPECT_ANY_THROW(testFunction->Evaluate(-1, at111_.data(), 0));
  }

  TEST_F(TranslatedFunctionTest, ShouldThrowWhenEvaluatedComponentTooLarge)
  {
    Teuchos::RCP<DRT::UTILS::Function> testFunction =
        Teuchos::rcp(new DRT::UTILS::TranslatedFunction(linearTime_, scalarField_));
    EXPECT_ANY_THROW(testFunction->Evaluate(2, at111_.data(), 0));
  }

  TEST_F(TranslatedFunctionTest, ShouldThrowWhenOriginHasOnlyOneComponent)
  {
    EXPECT_ANY_THROW(new DRT::UTILS::TranslatedFunction(scalarField_, constant_));
  }
}  // namespace