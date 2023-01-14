/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the muscle utils
\level 3
*/
/*----------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "mat_muscle_utils.H"
#include "unittests_assertions.h"

namespace
{
  TEST(MuscleUtilsTest, TestEvaluateLambert)
  {
    const double xi = 1.39234;
    double W0 = 0.89437;
    const double tol = 1e-15;
    const double maxiter = 100;

    const double ref_W = 0.69492997657856426;

    MAT::UTILS::MUSCLE::EvaluateLambert(xi, W0, tol, maxiter);

    EXPECT_NEAR(W0, ref_W, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluateForceStretchDependencyEhret)
  {
    const double l_smaller_lmin = 0.5;
    const double l_greater_lmin_smaller_lopt = 0.7;
    const double l_greater_lopt = 1.4;

    const double lmin = 0.6;
    const double lopt = 1.2;

    double ref_fxi_l_smaller_lmin = 0.0;
    double ref_fxi_l_equal_lmin = 0.0;
    double ref_fxi_l_greater_lmin_smaller_lopt = 0.27099677511525644;
    double ref_fxi_l_equal_lopt = 1.0;
    double ref_fxi_l_greater_lopt = 0.90374610400726718;

    auto test_fxi_l_smaller_lmin =
        MAT::UTILS::MUSCLE::EvaluateForceStretchDependencyEhret(l_smaller_lmin, lmin, lopt);
    auto test_fxi_l_equal_lmin =
        MAT::UTILS::MUSCLE::EvaluateForceStretchDependencyEhret(lmin, lmin, lopt);
    auto test_fxi_l_greater_lmin_smaller_lopt =
        MAT::UTILS::MUSCLE::EvaluateForceStretchDependencyEhret(
            l_greater_lmin_smaller_lopt, lmin, lopt);
    auto test_fxi_l_equal_lopt =
        MAT::UTILS::MUSCLE::EvaluateForceStretchDependencyEhret(lopt, lmin, lopt);
    auto test_fxi_l_greater_lopt =
        MAT::UTILS::MUSCLE::EvaluateForceStretchDependencyEhret(l_greater_lopt, lmin, lopt);

    EXPECT_NEAR(test_fxi_l_smaller_lmin, ref_fxi_l_smaller_lmin, 1.0e-10);
    EXPECT_NEAR(test_fxi_l_equal_lmin, ref_fxi_l_equal_lmin, 1.0e-10);
    EXPECT_NEAR(test_fxi_l_greater_lmin_smaller_lopt, ref_fxi_l_greater_lmin_smaller_lopt, 1.0e-10);
    EXPECT_NEAR(test_fxi_l_equal_lopt, ref_fxi_l_equal_lopt, 1.0e-10);
    EXPECT_NEAR(test_fxi_l_greater_lopt, ref_fxi_l_greater_lopt, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluateDerivativeForceStretchDependencyEhret)
  {
    const double l_smaller_lmin = 0.5;
    const double l_greater_lmin_smaller_lopt = 0.7;
    const double l_greater_lopt = 1.4;

    const double lmin = 0.6;
    const double lopt = 1.2;

    double ref_dfxi_l_smaller_lmin = 0.0;
    double ref_dfxi_l_equal_lmin = 0.0;
    double ref_dfxi_l_greater_lmin_smaller_lopt = 2.6346908691761053;
    double ref_dfxi_l_equal_lopt = 0.0;
    double ref_dfxi_l_greater_lopt = -0.87864204556262071;

    auto test_dfxi_l_smaller_lmin =
        MAT::UTILS::MUSCLE::EvaluateDerivativeForceStretchDependencyEhret(
            l_smaller_lmin, lmin, lopt);
    auto test_dfxi_l_equal_lmin =
        MAT::UTILS::MUSCLE::EvaluateDerivativeForceStretchDependencyEhret(lmin, lmin, lopt);
    auto test_dfxi_l_greater_lmin_smaller_lopt =
        MAT::UTILS::MUSCLE::EvaluateDerivativeForceStretchDependencyEhret(
            l_greater_lmin_smaller_lopt, lmin, lopt);
    auto test_dfxi_l_equal_lopt =
        MAT::UTILS::MUSCLE::EvaluateDerivativeForceStretchDependencyEhret(lopt, lmin, lopt);
    auto test_dfxi_l_greater_lopt =
        MAT::UTILS::MUSCLE::EvaluateDerivativeForceStretchDependencyEhret(
            l_greater_lopt, lmin, lopt);

    EXPECT_NEAR(test_dfxi_l_smaller_lmin, ref_dfxi_l_smaller_lmin, 1.0e-10);
    EXPECT_NEAR(test_dfxi_l_equal_lmin, ref_dfxi_l_equal_lmin, 1.0e-10);
    EXPECT_NEAR(
        test_dfxi_l_greater_lmin_smaller_lopt, ref_dfxi_l_greater_lmin_smaller_lopt, 1.0e-10);
    EXPECT_NEAR(test_dfxi_l_equal_lopt, ref_dfxi_l_equal_lopt, 1.0e-10);
    EXPECT_NEAR(test_dfxi_l_greater_lopt, ref_dfxi_l_greater_lopt, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluateForceVelocityDependencyBoel)
  {
    const double l_case_dotl_greater_zero = 1.2;
    const double l_case_dotl_smaller_zero = 0.5;

    const double lold = 0.7;
    const double timestepsize = 0.1;
    const double dotlmin = 0.5;
    const double de = 0.7;
    const double dc = 0.5;
    const double ke = 0.25;
    const double kc = 0.35;

    double ref_case_dotl_greater_zero = -59.9;
    double ref_case_dotl_smaller_zero = -5.75;

    auto test_case_dotl_greater_zero = MAT::UTILS::MUSCLE::EvaluateForceVelocityDependencyBoel(
        l_case_dotl_greater_zero, lold, timestepsize, dotlmin, de, dc, ke, kc);
    auto test_case_dotl_smaller_zero = MAT::UTILS::MUSCLE::EvaluateForceVelocityDependencyBoel(
        l_case_dotl_smaller_zero, lold, timestepsize, dotlmin, de, dc, ke, kc);

    EXPECT_NEAR(test_case_dotl_greater_zero, ref_case_dotl_greater_zero, 1.0e-10);
    EXPECT_NEAR(test_case_dotl_smaller_zero, ref_case_dotl_smaller_zero, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluateDerivativeForceVelocityDependencyBoel)
  {
    const double l_case_dotl_greater_zero = 1.2;
    const double l_case_dotl_smaller_zero = 0.5;

    const double lold = 0.7;
    const double timestepsize = 0.1;
    const double dotlmin = 0.5;
    const double de = 0.7;
    const double dc = 0.5;
    const double ke = 0.25;
    const double kc = 0.35;

    double ref_case_dotl_greater_zero = -974.4;
    double ref_case_dotl_smaller_zero = -84.375;

    auto test_case_dotl_greater_zero =
        MAT::UTILS::MUSCLE::EvaluateDerivativeForceVelocityDependencyBoel(
            l_case_dotl_greater_zero, lold, timestepsize, dotlmin, de, dc, ke, kc);
    auto test_case_dotl_smaller_zero =
        MAT::UTILS::MUSCLE::EvaluateDerivativeForceVelocityDependencyBoel(
            l_case_dotl_smaller_zero, lold, timestepsize, dotlmin, de, dc, ke, kc);

    EXPECT_NEAR(test_case_dotl_greater_zero, ref_case_dotl_greater_zero, 1.0e-10);
    EXPECT_NEAR(test_case_dotl_smaller_zero, ref_case_dotl_smaller_zero, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluateTimeDependentActiveStressEhret)
  {
    const double Na = 48.5435;
    const int muTypesNum = 3;
    const std::vector<double> rho = {0.4, 0.35, 0.25};
    const std::vector<double> I = {0.02, 0.04, 0.06};
    const std::vector<double> F = {0.3, 0.2, 0.1};
    const std::vector<double> T = {0.01, 0.02, 0.03};

    const int actIntervalsNum = 2;
    const std::vector<double> actTimes = {0.0, 0.1, 0.5};
    const std::vector<double> actValues = {0.0, 1.0};

    const double t_0 = 0.0;
    const double t_smaller_tact = 0.05;
    const double t_greater_tact = 0.4;

    double ref_act_stress_t_0 = 0.0;
    double ref_act_stress_t_smaller_tact = 0.0;
    double ref_act_stress_t_greater_tact = 5.3471385137375966;

    auto test_act_stress_t_0 = MAT::UTILS::MUSCLE::EvaluateTimeDependentActiveStressEhret(
        Na, muTypesNum, rho, I, F, T, actIntervalsNum, actTimes, actValues, t_0);
    auto test_act_stress_t_smaller_tact =
        MAT::UTILS::MUSCLE::EvaluateTimeDependentActiveStressEhret(
            Na, muTypesNum, rho, I, F, T, actIntervalsNum, actTimes, actValues, t_smaller_tact);
    auto test_act_stress_t_greater_tact =
        MAT::UTILS::MUSCLE::EvaluateTimeDependentActiveStressEhret(
            Na, muTypesNum, rho, I, F, T, actIntervalsNum, actTimes, actValues, t_greater_tact);

    EXPECT_NEAR(test_act_stress_t_0, ref_act_stress_t_0, 1.0e-10);
    EXPECT_NEAR(test_act_stress_t_smaller_tact, ref_act_stress_t_smaller_tact, 1.0e-10);
    EXPECT_NEAR(test_act_stress_t_greater_tact, ref_act_stress_t_greater_tact, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluateActiveForceStretchDependencyBlemker)
  {
    const double l_smaller_lopt = 0.5;
    const double l_greater_lopt = 1.4;

    const double lopt = 1.2;

    double ref_fxi_l_smaller_lopt = 0.0025;
    double ref_fxi_l_equal_lopt = 1.0;
    double ref_fxi_l_greater_lopt = 0.8888888888888888;

    auto test_fxi_l_smaller_lopt =
        MAT::UTILS::MUSCLE::EvaluateActiveForceStretchDependencyBlemker(l_smaller_lopt, lopt);
    auto test_fxi_l_equal_lopt =
        MAT::UTILS::MUSCLE::EvaluateActiveForceStretchDependencyBlemker(lopt, lopt);
    auto test_fxi_l_greater_lopt =
        MAT::UTILS::MUSCLE::EvaluateActiveForceStretchDependencyBlemker(l_greater_lopt, lopt);

    EXPECT_NEAR(test_fxi_l_smaller_lopt, ref_fxi_l_smaller_lopt, 1.0e-10);
    EXPECT_NEAR(test_fxi_l_equal_lopt, ref_fxi_l_equal_lopt, 1.0e-10);
    EXPECT_NEAR(test_fxi_l_greater_lopt, ref_fxi_l_greater_lopt, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluateDerivativeActiveForceStretchDependencyBlemker)
  {
    const double l_smaller_lopt = 0.5;
    const double l_greater_lopt = 1.4;

    const double lopt = 1.2;

    double ref_dfxi_l_smaller_lopt = 0.25;
    double ref_dfxi_l_equal_lopt = 0.0;
    double ref_dfxi_l_greater_lopt = -1.1111111111111116;

    auto test_dfxi_l_smaller_lopt =
        MAT::UTILS::MUSCLE::EvaluateDerivativeActiveForceStretchDependencyBlemker(
            l_smaller_lopt, lopt);
    auto test_dfxi_l_equal_lopt =
        MAT::UTILS::MUSCLE::EvaluateDerivativeActiveForceStretchDependencyBlemker(lopt, lopt);
    auto test_dfxi_l_greater_lopt =
        MAT::UTILS::MUSCLE::EvaluateDerivativeActiveForceStretchDependencyBlemker(
            l_greater_lopt, lopt);

    EXPECT_NEAR(test_dfxi_l_smaller_lopt, ref_dfxi_l_smaller_lopt, 1.0e-10);
    EXPECT_NEAR(test_dfxi_l_equal_lopt, ref_dfxi_l_equal_lopt, 1.0e-10);
    EXPECT_NEAR(test_dfxi_l_greater_lopt, ref_dfxi_l_greater_lopt, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluatePassiveForceStretchDependencyBlemker)
  {
    const double l_smaller_lopt = 0.8;
    const double l_greater_lopt_smaller_lstar = 1.3;
    const double l_greater_lstar = 1.5;

    const double lopt = 0.9;
    const double lstar = 1.4;

    const double P1 = 3.4832;
    const double P2 = 1.3234;

    double ref_fxi_l_smaller_lopt = 0.0;
    double ref_fxi_l_equal_lopt = 0.0;
    double ref_fxi_l_greater_lopt_smaller_lstar = 2.7890126634661554;
    double ref_fxi_l_equal_lstar = 3.7825653114188338;
    double ref_fxi_l_greater_lstar = 4.9696657821758325;

    auto test_fxi_l_smaller_lopt = MAT::UTILS::MUSCLE::EvaluatePassiveForceStretchDependencyBlemker(
        l_smaller_lopt, lopt, lstar, P1, P2);
    auto test_fxi_l_equal_lopt =
        MAT::UTILS::MUSCLE::EvaluatePassiveForceStretchDependencyBlemker(lopt, lopt, lstar, P1, P2);
    auto test_fxi_l_greater_lopt_smaller_lstar =
        MAT::UTILS::MUSCLE::EvaluatePassiveForceStretchDependencyBlemker(
            l_greater_lopt_smaller_lstar, lopt, lstar, P1, P2);
    auto test_fxi_l_equal_lstar = MAT::UTILS::MUSCLE::EvaluatePassiveForceStretchDependencyBlemker(
        lstar, lopt, lstar, P1, P2);
    auto test_fxi_l_greater_lstar =
        MAT::UTILS::MUSCLE::EvaluatePassiveForceStretchDependencyBlemker(
            l_greater_lstar, lopt, lstar, P1, P2);

    EXPECT_NEAR(test_fxi_l_smaller_lopt, ref_fxi_l_smaller_lopt, 1.0e-10);
    EXPECT_NEAR(test_fxi_l_equal_lopt, ref_fxi_l_equal_lopt, 1.0e-10);
    EXPECT_NEAR(
        test_fxi_l_greater_lopt_smaller_lstar, ref_fxi_l_greater_lopt_smaller_lstar, 1.0e-10);
    EXPECT_NEAR(test_fxi_l_equal_lstar, ref_fxi_l_equal_lstar, 1.0e-10);
    EXPECT_NEAR(test_fxi_l_greater_lstar, ref_fxi_l_greater_lstar, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluateDerivativePassiveForceStretchDependencyBlemker)
  {
    const double l_smaller_lopt = 0.8;
    const double l_greater_lopt_smaller_lstar = 1.3;
    const double l_greater_lstar = 1.5;

    const double lopt = 0.9;
    const double lstar = 1.4;

    const double P1 = 3.4832;
    const double P2 = 1.3234;

    double ref_dfxi_l_smaller_lopt = 0.0;
    double ref_dfxi_l_equal_lopt = 0.0;
    double ref_dfxi_l_greater_lopt_smaller_lstar = 9.2229402653678996;
    double ref_dfxi_l_equal_lstar = 11.871004707569982;
    double ref_dfxi_l_greater_lstar = 11.871004707569982;

    auto test_dfxi_l_smaller_lopt =
        MAT::UTILS::MUSCLE::EvaluateDerivativePassiveForceStretchDependencyBlemker(
            l_smaller_lopt, lopt, lstar, P1, P2);
    auto test_dfxi_l_equal_lopt =
        MAT::UTILS::MUSCLE::EvaluateDerivativePassiveForceStretchDependencyBlemker(
            lopt, lopt, lstar, P1, P2);
    auto test_dfxi_l_greater_lopt_smaller_lstar =
        MAT::UTILS::MUSCLE::EvaluateDerivativePassiveForceStretchDependencyBlemker(
            l_greater_lopt_smaller_lstar, lopt, lstar, P1, P2);
    auto test_dfxi_l_equal_lstar =
        MAT::UTILS::MUSCLE::EvaluateDerivativePassiveForceStretchDependencyBlemker(
            lstar, lopt, lstar, P1, P2);
    auto test_dfxi_l_greater_lstar =
        MAT::UTILS::MUSCLE::EvaluateDerivativePassiveForceStretchDependencyBlemker(
            l_greater_lstar, lopt, lstar, P1, P2);

    EXPECT_NEAR(test_dfxi_l_smaller_lopt, ref_dfxi_l_smaller_lopt, 1.0e-10);
    EXPECT_NEAR(test_dfxi_l_equal_lopt, ref_dfxi_l_equal_lopt, 1.0e-10);
    EXPECT_NEAR(
        test_dfxi_l_greater_lopt_smaller_lstar, ref_dfxi_l_greater_lopt_smaller_lstar, 1.0e-10);
    EXPECT_NEAR(test_dfxi_l_equal_lstar, ref_dfxi_l_equal_lstar, 1.0e-10);
    EXPECT_NEAR(test_dfxi_l_greater_lstar, ref_dfxi_l_greater_lstar, 1.0e-10);
  }

  TEST(MuscleUtilsTest, TestEvaluateTimeDependentActiveStressTanh)
  {
    const double sigma_max = 56.789;
    const double alpha = 1.1;
    const double beta = 34.567;
    const double t_act_start = 0.1;

    const double t_0 = 0.0;
    const double t_smaller_tact = 0.05;
    const double t_greater_tact = 0.4;

    double ref_act_stress_t_0 = 0.0;
    double ref_act_stress_t_smaller_tact = 0.0;
    double ref_act_stress_t_greater_tact = 62.467899877162083;

    auto test_act_stress_t_0 = MAT::UTILS::MUSCLE::EvaluateTimeDependentActiveStressTanh(
        sigma_max, alpha, beta, t_act_start, t_0);
    auto test_act_stress_t_smaller_tact = MAT::UTILS::MUSCLE::EvaluateTimeDependentActiveStressTanh(
        sigma_max, alpha, beta, t_act_start, t_smaller_tact);
    auto test_act_stress_t_greater_tact = MAT::UTILS::MUSCLE::EvaluateTimeDependentActiveStressTanh(
        sigma_max, alpha, beta, t_act_start, t_greater_tact);

    EXPECT_NEAR(test_act_stress_t_0, ref_act_stress_t_0, 1.0e-10);
    EXPECT_NEAR(test_act_stress_t_smaller_tact, ref_act_stress_t_smaller_tact, 1.0e-10);
    EXPECT_NEAR(test_act_stress_t_greater_tact, ref_act_stress_t_greater_tact, 1.0e-10);
  }
}  // namespace