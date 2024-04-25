/*----------------------------------------------------------------------*/
/*! \file

\brief unittests for LungOxygenExchangeLaw class in poromultiphase_scatra-framework

\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_poromultiphase_scatra_function.hpp"

#include <memory>
#include <string>
#include <vector>


namespace
{
  using namespace FourC;

  class LungOxygenExchangeLawTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      // function parameters
      const std::vector<std::pair<std::string, double>> func_params = {{"rho_oxy", 1.429e-9},
          {"DiffAdVTLC", 5.36}, {"alpha_oxy", 2.1e-4}, {"rho_air", 1.0e-9}, {"rho_bl", 1.03e-6},
          {"n", 3}, {"P_oB50", 3.6}, {"NC_Hb", 0.25}, {"P_atmospheric", 101.3},
          {"volfrac_blood_ref", 0.1}};

      // construct LungOxygenExchangeLaw
      LungOxygenExchangeLaw_ =
          std::make_unique<POROMULTIPHASESCATRA::LungOxygenExchangeLaw<3>>(func_params);
    }

    std::unique_ptr<POROMULTIPHASESCATRA::LungOxygenExchangeLaw<3>> LungOxygenExchangeLaw_;
  };

  class LungCarbonDioxideExchangeLawTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      // function parameters
      const std::vector<std::pair<std::string, double>> func_params = {{"rho_CO2", 1.98e-9},
          {"DiffsolAdVTLC", 4.5192e-3}, {"pH", 7.352}, {"rho_air", 1.0e-9}, {"rho_bl", 1.03e-6},
          {"rho_oxy", 1.429e-9}, {"n", 3}, {"P_oB50", 3.6}, {"C_Hb", 18.2}, {"NC_Hb", 0.25},
          {"alpha_oxy", 2.1e-4}, {"P_atmospheric", 101.3}, {"ScalingFormmHg", 133.3e-3},
          {"volfrac_blood_ref", 0.1}};

      // construct LungCarbonDioxideExchangeLaw
      LungCarbonDioxideExchangeLaw_ =
          std::make_unique<POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<3>>(func_params);
    }

    std::unique_ptr<POROMULTIPHASESCATRA::LungCarbonDioxideExchangeLaw<3>>
        LungCarbonDioxideExchangeLaw_;
  };

  TEST_F(LungOxygenExchangeLawTest, TestEvaluateAndEvaluateDerivativeZeroOxygenatedBlood)
  {
    // input arguments
    const int component = 0;
    const std::vector<std::pair<std::string, double>> variables = {{"phi1", 0.19}, {"phi2", 0.0}};
    const std::vector<std::pair<std::string, double>> constants = {
        {"p1", 0.005}, {"S1", 0.005}, {"porosity", 0.2}, {"VF1", 0.3}};

    // test Evaluate
    EXPECT_NEAR(LungOxygenExchangeLaw_->Evaluate(variables, constants, component),
        6.4996477560000005e-11, 1e-14);

    // test EvaluateDerivative wrt phi1
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(variables, constants, component)[0],
        3.42086724e-10, 1e-14);

    // test EvaluateDerivative wrt phi2
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(variables, constants, component)[1],
        -1.0016939520000001e-07, 1e-14);
  }

  TEST_F(LungOxygenExchangeLawTest, TestEvaluateAndEvaluateDerivativeHalfOxygenatedBlood)
  {
    // input arguments
    const int component = 0;
    const std::vector<std::pair<std::string, double>> variables = {
        {"phi1", 0.19}, {"phi2", 1.5e-4}};
    const std::vector<std::pair<std::string, double>> constants = {
        {"p1", 0.005}, {"S1", 0.005}, {"porosity", 0.2}, {"VF1", 0.2}};

    // test Evaluate
    EXPECT_NEAR(LungOxygenExchangeLaw_->Evaluate(variables, constants, component),
        3.2792446508136922e-11, 1e-14);

    // test EvaluateDerivative wrt phi1
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(variables, constants, component)[0],
        2.2805781600000004e-10, 1e-14);

    // test EvaluateDerivative wrt phi2
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(variables, constants, component)[1],
        -4.1174493452973438e-08, 1e-14);
  }

  TEST_F(LungOxygenExchangeLawTest, TestEvaluateAndEvaluateDerivativeNearlyFullyOxygenatedBlood)
  {
    // input arguments
    const int component = 0;
    const std::vector<std::pair<std::string, double>> variables = {
        {"phi1", 0.19}, {"phi2", 3.0e-4}};
    const std::vector<std::pair<std::string, double>> constants = {
        {"p1", 0.005}, {"S1", 0.005}, {"porosity", 0.2}, {"VF1", 0.03}};

    // test Evaluate
    EXPECT_NEAR(LungOxygenExchangeLaw_->Evaluate(variables, constants, component),
        3.3233325598255571e-12, 1e-14);

    // test EvaluateDerivative wrt phi1
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(variables, constants, component)[0],
        3.4208672400000007e-11, 1e-14);

    // test EvaluateDerivative wrt phi2
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(variables, constants, component)[1],
        -2.4885189707947446e-08, 1e-14);
  }

  TEST_F(LungCarbonDioxideExchangeLawTest,
      TestEvaluateAndEvaluateDerivativeNearlyFullyCarbonDioxideEnrichedBlood)
  {
    // input arguments
    const int component = 0;
    const std::vector<std::pair<std::string, double>> variables = {
        {"phi1", 0.19}, {"phi2", 1.88e-4}, {"phi3", 0.06}, {"phi4", 0.0894}};
    const std::vector<std::pair<std::string, double>> constants = {
        {"p1", 0.005}, {"S1", 0.0}, {"porosity", 0.0}, {"VF1", 0.4}};

    // test Evaluate
    EXPECT_NEAR(LungCarbonDioxideExchangeLaw_->Evaluate(variables, constants, component),
        1.0687549509499465e-10, 1e-14);

    // test EvaluateDerivative wrt phi1
    EXPECT_NEAR(
        LungCarbonDioxideExchangeLaw_->EvaluateDerivative(variables, constants, component)[2],
        -1.8312702239999997e-09, 1e-14);

    // test EvaluateDerivative wrt phi2
    EXPECT_NEAR(
        LungCarbonDioxideExchangeLaw_->EvaluateDerivative(variables, constants, component)[3],
        2.4245157554249963e-09, 1e-14);
  }

  TEST_F(LungCarbonDioxideExchangeLawTest, TestEvaluateAndEvaluateDerivativeCarbonDioxidePoorBlood)
  {
    // input arguments
    const int component = 0;
    const std::vector<std::pair<std::string, double>> variables = {
        {"phi1", 0.19}, {"phi2", 1.88e-4}, {"phi3", 0.06}, {"phi4", 0.08760}};
    const std::vector<std::pair<std::string, double>> constants = {
        {"p1", 0.005}, {"S1", 0.0}, {"porosity", 0.0}, {"VF1", 0.4}};

    // test Evaluate
    EXPECT_NEAR(LungCarbonDioxideExchangeLaw_->Evaluate(variables, constants, component),
        1.0251136673522964e-10, 1e-14);

    // test EvaluateDerivative wrt phi1
    EXPECT_NEAR(
        LungCarbonDioxideExchangeLaw_->EvaluateDerivative(variables, constants, component)[2],
        -1.8312702239999997e-09, 1e-14);

    // test EvaluateDerivative wrt phi2
    EXPECT_NEAR(
        LungCarbonDioxideExchangeLaw_->EvaluateDerivative(variables, constants, component)[3],
        2.4245157554249963e-09, 1e-14);
  }

}  // namespace