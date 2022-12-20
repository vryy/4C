/*----------------------------------------------------------------------*/
/*! \file

\brief unittests for LungOxygenExchangeLaw class in poromultiphase_scatra-framework

\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "poromultiphase_scatra_function.H"

#include <vector>
#include <string>
#include <memory>


namespace
{
  class LungOxygenExchangeLawTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      // function parameters
      const std::vector<std::pair<std::string, double>> func_params = {{"rho_oxy", 1.429e-9},
          {"DiffAdVTLC", 5.36}, {"alpha_oxy", 2.1e-4}, {"rho_air", 1.0e-9}, {"rho_bl", 1.03e-6},
          {"n", 3}, {"P_oB50", 3.6}, {"NC_Hb", 0.25}, {"P_atmospheric", 101.3}};

      // construct LungOxygenExchangeLaw
      LungOxygenExchangeLaw_ = std::unique_ptr<POROMULTIPHASESCATRA::LungOxygenExchangeLaw<3>>(
          new POROMULTIPHASESCATRA::LungOxygenExchangeLaw<3>(func_params));
    }

    std::unique_ptr<POROMULTIPHASESCATRA::LungOxygenExchangeLaw<3>> LungOxygenExchangeLaw_;
  };

  TEST_F(LungOxygenExchangeLawTest, TestEvaluateAndEvaluateDerivativeZeroOxygenatedBlood)
  {
    // input arguments
    const int index = 0;
    const std::vector<std::pair<std::string, double>> variables = {{"phi1", 0.19}, {"phi2", 0.0}};
    const std::vector<std::pair<std::string, double>> constants = {{"p1", 0.005}};

    // test Evaluate
    EXPECT_NEAR(
        LungOxygenExchangeLaw_->Evaluate(index, variables, constants), 2.166549252e-11, 1e-14);

    // test EvaluateDerivative wrt phi1
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(index, variables, constants)[0],
        1.14028908e-10, 1e-14);

    // test EvaluateDerivative wrt phi2
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(index, variables, constants)[1],
        -3.33897984e-08, 1e-14);
  }

  TEST_F(LungOxygenExchangeLawTest, TestEvaluateAndEvaluateDerivativeHalfOxygenatedBlood)
  {
    // input arguments
    const int index = 0;
    const std::vector<std::pair<std::string, double>> variables = {
        {"phi1", 0.19}, {"phi2", 1.5e-4}};
    const std::vector<std::pair<std::string, double>> constants = {{"p1", 0.005}};

    // test Evaluate
    EXPECT_NEAR(
        LungOxygenExchangeLaw_->Evaluate(index, variables, constants), 1.639622325407e-11, 1e-14);

    // test EvaluateDerivative wrt phi1
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(index, variables, constants)[0],
        1.14028908e-10, 1e-14);

    // test EvaluateDerivative wrt phi2
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(index, variables, constants)[1],
        -2.058724672649e-08, 1e-14);
  }

  TEST_F(LungOxygenExchangeLawTest, TestEvaluateAndEvaluateDerivativeNearlyFullyOxygenatedBlood)
  {
    // input arguments
    const int index = 0;
    const std::vector<std::pair<std::string, double>> variables = {
        {"phi1", 0.19}, {"phi2", 3.0e-4}};
    const std::vector<std::pair<std::string, double>> constants = {{"p1", 0.005}};

    // test Evaluate
    EXPECT_NEAR(
        LungOxygenExchangeLaw_->Evaluate(index, variables, constants), 1.10777752e-11, 1e-14);

    // test EvaluateDerivative wrt phi1
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(index, variables, constants)[0],
        1.14028908e-10, 1e-14);

    // test EvaluateDerivative wrt phi2
    EXPECT_NEAR(LungOxygenExchangeLaw_->EvaluateDerivative(index, variables, constants)[1],
        -8.295063236e-08, 1e-14);
  }

}  // namespace