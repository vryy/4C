/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the full constrained mixture fiber
\level 3
*/
/*----------------------------------------------------------------------*/

#include <memory>
#include <gtest/gtest.h>

#include <Sacado.hpp>

#include "mat_par_material.H"
#include "mixture_growth_evolution_linear_cauchy_poisson_turnover.H"
#include "mixture_constituent_remodelfiber_material_exponential.H"
#include "mixture_full_constrained_mixture_fiber.H"

#include "linalg_fixedsizematrix.H"

#include "mixture_remodelfiber-internal.H"

namespace
{
  using FADdouble = Sacado::Fad::DFad<double>;
  class FullConstrainedMixtureFiberTest : public ::testing::Test
  {
   protected:
    template <typename Number>
    const MIXTURE::PAR::RemodelFiberMaterialExponential<Number>* CreateMaterial()
    {
      const auto container = Teuchos::rcp(new MAT::PAR::Material());
      container->Add("K1", 1.3);
      container->Add("K2", 1.3);
      container->Add("COMPRESSION", true);
      static auto params =
          std::make_shared<MIXTURE::PAR::RemodelFiberMaterialExponential<Number>>(container);

      return params.get();
    }

    template <typename Number>
    MIXTURE::FullConstrainedMixtureFiber<Number> GenerateFiber(const double decay_time = 12.0,
        const double growth_constant = 0.1, const Number lambda_pre = 1.1,
        const bool adaptive = false, const bool growth_enabled = true)
    {
      return {std::make_shared<MIXTURE::RemodelFiberMaterialExponential<Number>>(
                  CreateMaterial<Number>()),
          {growth_constant, decay_time}, lambda_pre, adaptive, growth_enabled};
    }

    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, double> GenerateRemodelFiber(
        const double decay_time = 12.0, const double growth_constant = 0.1,
        const double lambda_pre = 1.1)
    {
      return {std::make_unique<const MIXTURE::RemodelFiberMaterialExponential<double>>(
                  CreateMaterial<double>()),
          {growth_constant, decay_time}, lambda_pre};
    }
  };

  TEST_F(FullConstrainedMixtureFiberTest, TestThrowsIfHistoryIsEmpty)
  {
#ifndef DEBUG
    GTEST_SKIP() << "Skip debug assertion tests in release mode.";
#endif
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>();

    const double lambda_f = 1.0;

    ASSERT_ANY_THROW(fiber.RecomputeState(lambda_f, 0.0));
  }

  TEST_F(FullConstrainedMixtureFiberTest, TestDoesNotThrowIfHistoryIsEmptyButGrowthIsDisabled)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>(12.0, 0.1, 1.1, true, false);

    const double lambda_f = 1.0;

    fiber.RecomputeState(lambda_f, 0.0);
  }

  TEST_F(FullConstrainedMixtureFiberTest, TestEvaluateDGrowthEvolutionEquationDtDGrowth)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>();
    const double lambda_f = 1.0;
    fiber.ReinitializeHistory(lambda_f, 0.0);

    fiber.RecomputeState(lambda_f, 0.0);

    ASSERT_NEAR(fiber.computed_sigma_, fiber.sig_h_, 1e-18);
  }

  TEST_F(FullConstrainedMixtureFiberTest, TestFiberMaintenanceAfterOneStep)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>();
    const double lambda_f = 1.0;
    fiber.ReinitializeHistory(lambda_f, 0.0);

    fiber.RecomputeState(lambda_f, 0.01);

    EXPECT_NEAR(fiber.computed_growth_scalar_, 1.0, 1e-9);
    EXPECT_NEAR(fiber.computed_sigma_, fiber.sig_h_, 1e-9);
  }

  TEST_F(FullConstrainedMixtureFiberTest, TestFiberMaintenanceTwoSteps)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>();
    const double lambda_f = 1.0;
    fiber.ReinitializeHistory(lambda_f, 0.0);

    for (unsigned timestep = 1; timestep <= 2; ++timestep)
    {
      fiber.RecomputeState(lambda_f, timestep * 0.01);
      fiber.Update();
    }

    fiber.RecomputeState(lambda_f, 0.1);

    EXPECT_NEAR(fiber.computed_growth_scalar_, 1.0, 1e-9);
    EXPECT_NEAR(fiber.computed_sigma_, fiber.sig_h_, 1e-9);
  }

  TEST_F(FullConstrainedMixtureFiberTest, TestFiberMaintenanceNSteps)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>();
    const double lambda_f = 1.0;
    fiber.ReinitializeHistory(lambda_f, 0.0);

    for (unsigned timestep = 1; timestep <= 3; ++timestep)
    {
      fiber.RecomputeState(lambda_f, timestep * 0.01);
      fiber.Update();
    }

    fiber.RecomputeState(lambda_f, 0.1);

    EXPECT_NEAR(fiber.computed_growth_scalar_, 1.0, 1e-9);
    EXPECT_NEAR(fiber.computed_sigma_, fiber.sig_h_, 1e-9);
  }

  TEST_F(FullConstrainedMixtureFiberTest, CheckDerivativeGrowthScalarIntegrand)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<FADdouble>();
    const MIXTURE::MassIncrement<FADdouble> mass_increment{
        1.01, FADdouble(2, 0, 1.12), FADdouble(2, 1, 1.12), 1.0};

    const double time = 1.4;
    const FADdouble growth_scalar_integrand = cm_fiber.GrowthScalarIntegrand(mass_increment, time);

    EXPECT_FLOAT_EQ(growth_scalar_integrand.dx(0),
        cm_fiber.DGrowthScalarIntegrandDGrowthScalar(mass_increment, time).val());

    EXPECT_FLOAT_EQ(growth_scalar_integrand.dx(1),
        cm_fiber.DGrowthScalarIntegrandDProductionRate(mass_increment, time).val());
  }

  TEST_F(FullConstrainedMixtureFiberTest, CheckDerivativeScaledCauchyStressIntegrand)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<FADdouble>();

    const double time = 1.4;
    const FADdouble lambda_f = 1.012;
    cm_fiber.ReinitializeHistory(lambda_f, time);
    cm_fiber.RecomputeState(lambda_f, time);
    const MIXTURE::MassIncrement<FADdouble> mass_increment{
        1.01, FADdouble(2, 0, 1.12), FADdouble(2, 1, 1.12), 1.0};

    const FADdouble scaled_cauchy_stress_integrand =
        cm_fiber.ScaledCauchyStressIntegrand(mass_increment, time, lambda_f);

    EXPECT_FLOAT_EQ(scaled_cauchy_stress_integrand.dx(0),
        cm_fiber.DScaledCauchyStressIntegrandDGrowthScalar(mass_increment, time, lambda_f).val());

    EXPECT_FLOAT_EQ(scaled_cauchy_stress_integrand.dx(1),
        cm_fiber.DScaledCauchyStressIntegrandDProductionRate(mass_increment, time, lambda_f).val());
  }

  TEST_F(FullConstrainedMixtureFiberTest, RemodelFiberComparisonCauchyStress)
  {
    // compare the results of the full constrained mixture model with the homogenized constrained
    // mixture model
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();
    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, double> remodel_fiber =
        GenerateRemodelFiber();

    const double lambda_f = 1.05;
    cm_fiber.ReinitializeHistory(lambda_f, 0.0);

    const double dt = 0.1;
    for (unsigned timestep = 1; timestep <= 1000; ++timestep)
    {
      const double time = timestep * dt;
      cm_fiber.RecomputeState(lambda_f, time);
      remodel_fiber.SetState(lambda_f, 1.0);

      remodel_fiber.IntegrateLocalEvolutionEquationsImplicit(dt);

      // compare
      ASSERT_NEAR(
          cm_fiber.sig_h_, remodel_fiber.EvaluateCurrentHomeostaticFiberCauchyStress(), 1e-2);
      ASSERT_NEAR(cm_fiber.computed_sigma_ / cm_fiber.sig_h_,
          remodel_fiber.EvaluateCurrentFiberCauchyStress() /
              remodel_fiber.EvaluateCurrentHomeostaticFiberCauchyStress(),
          1e-2);
      ASSERT_NEAR(
          cm_fiber.computed_growth_scalar_, remodel_fiber.EvaluateCurrentGrowthScalar(), 1e-2);

      cm_fiber.Update();
      remodel_fiber.Update();
    }
  }

  TEST_F(FullConstrainedMixtureFiberTest, RemodelFiberComparison2PKStress)
  {
    // compare the results of the full constrained mixture model with the homogenized constrained
    // mixture model
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();
    MIXTURE::FullConstrainedMixtureFiber cm_fiber_adaptive =
        GenerateFiber<double>(12.0, 0.1, 1.1, true);
    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, double> remodel_fiber =
        GenerateRemodelFiber();

    const double lambda_f = 1.05;
    cm_fiber.ReinitializeHistory(lambda_f, 0.0);
    cm_fiber_adaptive.ReinitializeHistory(lambda_f, 0.0);

    const double dt = 0.1;
    for (unsigned timestep = 1; timestep <= 1000; ++timestep)
    {
      const double time = timestep * dt;
      cm_fiber.RecomputeState(lambda_f, time);
      cm_fiber_adaptive.RecomputeState(lambda_f, time);

      remodel_fiber.SetState(lambda_f, 1.0);
      remodel_fiber.IntegrateLocalEvolutionEquationsImplicit(dt);

      // compare
      ASSERT_NEAR(cm_fiber.EvaluateCurrentSecondPKStress(),
          remodel_fiber.EvaluateCurrentFiberPK2Stress(), 1e-2);
      ASSERT_NEAR(cm_fiber_adaptive.EvaluateCurrentSecondPKStress(),
          remodel_fiber.EvaluateCurrentFiberPK2Stress(), 1e-2);

      cm_fiber.Update();
      cm_fiber_adaptive.Update();
      remodel_fiber.Update();
    }

    EXPECT_LE(cm_fiber_adaptive.history_[0].size(), 0.1 * cm_fiber.history_[0].size());
  }

  TEST_F(FullConstrainedMixtureFiberTest, LocalNewtonHasAnalyticalDerivativeInFirstTimestep)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<FADdouble>();

    const double lambda_f = 1.05;
    cm_fiber.ReinitializeHistory(lambda_f, 0.0);
    cm_fiber.RecomputeState(1.06, 0.01);

    // evaluate local newton
    const auto local_newton_evaluator = cm_fiber.GetLocalNewtonEvaluator();

    CORE::LINALG::Matrix<2, 1, FADdouble> x;
    x(0) = FADdouble(2, 0, 1.5);
    x(1) = FADdouble(2, 1, 0.453);

    const auto [residuum, derivative] = local_newton_evaluator(x);

    EXPECT_DOUBLE_EQ(residuum(0).dx(0), derivative(0, 0).val());
    EXPECT_DOUBLE_EQ(residuum(0).dx(1), derivative(0, 1).val());
    EXPECT_DOUBLE_EQ(residuum(1).dx(0), derivative(1, 0).val());
    EXPECT_DOUBLE_EQ(residuum(1).dx(1), derivative(1, 1).val());
  }

  TEST_F(FullConstrainedMixtureFiberTest, LocalNewtonHasAnalyticalDerivativeInFirstThreeTimesteps)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<FADdouble>();

    const double lambda_f = 1.05;
    cm_fiber.ReinitializeHistory(lambda_f, 0.0);
    const double dt = 0.1;
    for (unsigned timestep = 1; timestep <= 3; ++timestep)
    {
      const double time = timestep * dt;
      cm_fiber.RecomputeState(lambda_f, time);

      // evaluate local newton
      const auto local_newton_evaluator = cm_fiber.GetLocalNewtonEvaluator();

      CORE::LINALG::Matrix<2, 1, FADdouble> x;
      x(0) = FADdouble(2, 0, cm_fiber.computed_growth_scalar_.val());
      x(1) = FADdouble(2, 1, cm_fiber.computed_sigma_.val());

      const auto [residuum, derivative] = local_newton_evaluator(x);

      EXPECT_DOUBLE_EQ(residuum(0).dx(0), derivative(0, 0).val());
      EXPECT_DOUBLE_EQ(residuum(0).dx(1), derivative(0, 1).val());
      EXPECT_DOUBLE_EQ(residuum(1).dx(0), derivative(1, 0).val());
      EXPECT_DOUBLE_EQ(residuum(1).dx(1), derivative(1, 1).val());

      cm_fiber.Update();
    }
  }

  TEST_F(FullConstrainedMixtureFiberTest, DScaledCauchyStressIntegrandDLambdaFSq)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<FADdouble>();

    const FADdouble lambda_f = FADdouble(1, 0, 1.2);
    const double time = 1.0;
    cm_fiber.ReinitializeHistory(lambda_f.val(), 0.0);

    cm_fiber.RecomputeState(lambda_f, time);
    MIXTURE::MassIncrement<FADdouble> current_increment = {1.1, 1.4, 0.9, 1.22};

    FADdouble scaled_cauchy_stress_integrand =
        cm_fiber.ScaledCauchyStressIntegrand(current_increment, time, lambda_f);

    FADdouble dscaled_cauchy_stress_integrandDLambdaFSq =
        cm_fiber.DScaledCauchyStressIntegrandDLambdaFSq(current_increment, time, lambda_f);

    EXPECT_DOUBLE_EQ(dscaled_cauchy_stress_integrandDLambdaFSq.val() * 2 * lambda_f.val(),
        scaled_cauchy_stress_integrand.dx(0));
  }



  TEST_F(FullConstrainedMixtureFiberTest, DerivativeOfResiduumWithRespectToLambdaF)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<FADdouble>(999999.0, 0.1, 1.1);

    const double lambda_f_0 = 1.05;
    cm_fiber.ReinitializeHistory(lambda_f_0, 0.0);

    const double lambda_f = 1.06;
    cm_fiber.RecomputeState(FADdouble(1, 0, lambda_f), 0.01);

    // evaluate local newton
    const auto local_newton_evaluator = cm_fiber.GetLocalNewtonEvaluator();

    cm_fiber.computed_growth_scalar_ = cm_fiber.computed_growth_scalar_.val();
    cm_fiber.computed_sigma_ = cm_fiber.computed_sigma_.val();

    CORE::LINALG::Matrix<2, 1, FADdouble> x;
    x(0) = cm_fiber.computed_growth_scalar_;
    x(1) = cm_fiber.computed_sigma_;

    const auto [residuum, derivative] = local_newton_evaluator(x);

    EXPECT_DOUBLE_EQ(
        residuum(0).dx(0), cm_fiber.EvaluateDResiduumGrowthScalarDLambdaFSq().val() * 2 * lambda_f);
    EXPECT_DOUBLE_EQ(
        residuum(1).dx(0), cm_fiber.EvaluateDResiduumCauchyStressDLambdaFSq().val() * 2 * lambda_f);
  }

  TEST_F(FullConstrainedMixtureFiberTest,
      StressResponseLinearizationDuringHomeostasisIsAnalyticalSolution)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<FADdouble>();

    const double lambda_f_0 = 1.1;
    cm_fiber.ReinitializeHistory(lambda_f_0, 0.0);
    const double dt = 0.1;
    for (unsigned timestep = 1; timestep <= 3; ++timestep)
    {
      const double time = timestep * dt;
      const double lambda_f = lambda_f_0;
      cm_fiber.RecomputeState(FADdouble(1, 0, lambda_f), time);

      EXPECT_NEAR(cm_fiber.EvaluateDCurrentFiberPK2StressDLambdafsq().val() * 2 * lambda_f,
          cm_fiber.EvaluateCurrentSecondPKStress().dx(0), 1e-7);

      // remove automatic differentiation type from history data
      cm_fiber.computed_growth_scalar_ = cm_fiber.computed_growth_scalar_.val();
      cm_fiber.computed_sigma_ = cm_fiber.computed_sigma_.val();
      cm_fiber.computed_dgrowth_scalar_dlambda_f_sq_ =
          cm_fiber.computed_dgrowth_scalar_dlambda_f_sq_.val();
      cm_fiber.computed_dsigma_dlambda_f_sq_ = cm_fiber.computed_dsigma_dlambda_f_sq_.val();
      cm_fiber.current_state_.lambda_f = cm_fiber.current_state_.lambda_f.val();
      cm_fiber.Update();
    }
  }

  TEST_F(FullConstrainedMixtureFiberTest, StressResponseLinearizationIsAnalyticalSolution)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<FADdouble>();

    const double lambda_f_0 = 1.05;
    cm_fiber.ReinitializeHistory(lambda_f_0, 0.0);
    const double dt = 0.1;
    for (unsigned timestep = 1; timestep <= 3; ++timestep)
    {
      const double time = timestep * dt;
      const double lambda_f = lambda_f_0 + 0.001 * timestep;
      cm_fiber.RecomputeState(FADdouble(1, 0, lambda_f), time);

      EXPECT_NEAR(cm_fiber.EvaluateDCurrentFiberPK2StressDLambdafsq().val() * 2 * lambda_f,
          cm_fiber.EvaluateCurrentSecondPKStress().dx(0), 1e-8);

      // remove automatic differentiation type from history data
      cm_fiber.computed_growth_scalar_ = cm_fiber.computed_growth_scalar_.val();
      cm_fiber.computed_sigma_ = cm_fiber.computed_sigma_.val();
      cm_fiber.computed_dgrowth_scalar_dlambda_f_sq_ =
          cm_fiber.computed_dgrowth_scalar_dlambda_f_sq_.val();
      cm_fiber.computed_dsigma_dlambda_f_sq_ = cm_fiber.computed_dsigma_dlambda_f_sq_.val();
      cm_fiber.current_state_.lambda_f = cm_fiber.current_state_.lambda_f.val();
      cm_fiber.Update();
    }
  }

  TEST_F(FullConstrainedMixtureFiberTest,
      StressResponseLinearizationIsAnalyticalSolutionInInitialGrowthFreePeriod)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber =
        GenerateFiber<FADdouble>(12.0, 0.1, 1.1, true, false);

    double lambda_f = 1.2;
    fiber.RecomputeState(FADdouble(1, 0, lambda_f), 1.1);

    EXPECT_NEAR(fiber.EvaluateDCurrentFiberPK2StressDLambdafsq().val() * 2 * lambda_f,
        fiber.EvaluateCurrentSecondPKStress().dx(0), 1e-8);
  }

  TEST_F(FullConstrainedMixtureFiberTest, GrowthScalarLinearizationIsAnalyticalSolution)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<FADdouble>();

    const double lambda_f_0 = 1.05;
    cm_fiber.ReinitializeHistory(lambda_f_0, 0.0);
    const double dt = 0.1;
    for (unsigned timestep = 1; timestep <= 3; ++timestep)
    {
      const double time = timestep * dt;
      const double lambda_f = lambda_f_0 + 0.001 * timestep;
      cm_fiber.RecomputeState(FADdouble(1, 0, lambda_f), time);

      EXPECT_NEAR(cm_fiber.computed_dgrowth_scalar_dlambda_f_sq_.val() * 2 * lambda_f,
          cm_fiber.computed_growth_scalar_.dx(0), 1e-8);

      // remove automatic differentiation type from history data
      cm_fiber.computed_growth_scalar_ = cm_fiber.computed_growth_scalar_.val();
      cm_fiber.computed_sigma_ = cm_fiber.computed_sigma_.val();
      cm_fiber.computed_dgrowth_scalar_dlambda_f_sq_ =
          cm_fiber.computed_dgrowth_scalar_dlambda_f_sq_.val();
      cm_fiber.computed_dsigma_dlambda_f_sq_ = cm_fiber.computed_dsigma_dlambda_f_sq_.val();
      cm_fiber.current_state_.lambda_f = cm_fiber.current_state_.lambda_f.val();
      cm_fiber.Update();
    }
  }

  TEST_F(FullConstrainedMixtureFiberTest, TestReinitializeIdenticalHistory)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();

    cm_fiber.ReinitializeHistory(1.0, 0.0);

    ASSERT_EQ(cm_fiber.history_.size(), 1);
    ASSERT_EQ(cm_fiber.history_.back().size(), 1);

    cm_fiber.ReinitializeHistory(1.0, 0.0);
    ASSERT_EQ(cm_fiber.history_.size(), 1);
    ASSERT_EQ(cm_fiber.history_.back().size(), 1);
  }

  TEST_F(FullConstrainedMixtureFiberTest, TestReinitializeThrowsIfTimeIsNotIdentical)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();

    cm_fiber.ReinitializeHistory(1.0, 0.0);
    ASSERT_EQ(cm_fiber.history_.size(), 1);
    ASSERT_EQ(cm_fiber.history_.back().size(), 1);

    ASSERT_ANY_THROW(cm_fiber.ReinitializeHistory(1.0, 1.0));
  }

  TEST_F(FullConstrainedMixtureFiberTest, TestReinitializeAddNewBlockOverwritesIfOldIsJustTimestep)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();

    cm_fiber.ReinitializeHistory(1.0, 0.0);
    ASSERT_EQ(cm_fiber.history_.size(), 1);
    ASSERT_EQ(cm_fiber.history_.back().size(), 1);

    cm_fiber.ReinitializeHistory(1.2, 0.0);

    ASSERT_EQ(cm_fiber.history_.size(), 1);
    ASSERT_EQ(cm_fiber.history_.back().size(), 1);
  }

  TEST_F(FullConstrainedMixtureFiberTest, TestReinitializeAddNewBlock)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();

    cm_fiber.ReinitializeHistory(1.0, 0.0);
    ASSERT_EQ(cm_fiber.history_.size(), 1);
    ASSERT_EQ(cm_fiber.history_.back().size(), 1);

    cm_fiber.RecomputeState(1.1, 1.0);
    cm_fiber.Update();
    ASSERT_EQ(cm_fiber.history_.size(), 1);
    ASSERT_EQ(cm_fiber.history_.back().size(), 2);

    cm_fiber.ReinitializeHistory(1.2, 1.0);
    ASSERT_EQ(cm_fiber.history_.size(), 2);
    ASSERT_EQ(cm_fiber.history_[0].size(), 2);
    ASSERT_EQ(cm_fiber.history_[1].size(), 1);
  }

  TEST_F(FullConstrainedMixtureFiberTest, ComparisonWithMultipleIntervals)
  {
    // compare the results of the full constrained mixture model with the homogenized constrained
    // mixture model
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();
    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, double> remodel_fiber =
        GenerateRemodelFiber();

    const std::array lambda_f = {1.05, 1.1, 1.09999};
    double time = 0.0;
    unsigned timestep = 1;

    for (unsigned interval_id = 0; interval_id < lambda_f.size(); ++interval_id)
    {
      cm_fiber.ReinitializeHistory(lambda_f[interval_id], time);

      const double dt = 0.1;
      for (; timestep <= 1000 * (interval_id + 1); ++timestep)
      {
        time = timestep * dt;
        cm_fiber.RecomputeState(lambda_f[interval_id], time);
        remodel_fiber.SetState(lambda_f[interval_id], 1.0);

        remodel_fiber.IntegrateLocalEvolutionEquationsImplicit(dt);

        // compare
        ASSERT_NEAR(cm_fiber.computed_sigma_ / cm_fiber.sig_h_,
            remodel_fiber.EvaluateCurrentFiberCauchyStress() /
                remodel_fiber.EvaluateCurrentHomeostaticFiberCauchyStress(),
            1e-2);
        ASSERT_NEAR(
            cm_fiber.computed_growth_scalar_, remodel_fiber.EvaluateCurrentGrowthScalar(), 3e-2);

        cm_fiber.Update();
        remodel_fiber.Update();
      }
    }
  }


  TEST_F(FullConstrainedMixtureFiberTest, FiberWithJumpInTimeWithSingleInterval)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber_timeshift = GenerateFiber<double>();

    cm_fiber_timeshift.ReinitializeHistory(1.0, 0.0);

    for (int timestep = 1; timestep <= 4; ++timestep)
    {
      cm_fiber_timeshift.RecomputeState(1.1 + 0.001 * timestep, 0.1 * timestep);

      ASSERT_DOUBLE_EQ(
          cm_fiber_timeshift.computed_growth_scalar_, cm_fiber_timeshift.computed_growth_scalar_);

      cm_fiber_timeshift.Update();
    }

    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.reference_time_, 0.0);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[0].deposition_time, 0.0);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[1].deposition_time, 0.1);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[2].deposition_time, 0.2);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[3].deposition_time, 0.3);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[4].deposition_time, 0.4);

    cm_fiber_timeshift.AddTime(1.0);

    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.reference_time_, 1.0);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[0].deposition_time, 1.0);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[1].deposition_time, 1.1);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[2].deposition_time, 1.2);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[3].deposition_time, 1.3);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_.back()[4].deposition_time, 1.4);
  }

  TEST_F(FullConstrainedMixtureFiberTest, FiberWithJumpInTimeWithMultipleIntervals)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber_timeshift = GenerateFiber<double>();

    cm_fiber_timeshift.ReinitializeHistory(1.0, 0.0);

    for (int timestep = 1; timestep <= 3; ++timestep)
    {
      cm_fiber_timeshift.RecomputeState(1.1 + 0.001 * timestep, 0.1 * timestep);

      ASSERT_DOUBLE_EQ(
          cm_fiber_timeshift.computed_growth_scalar_, cm_fiber_timeshift.computed_growth_scalar_);

      cm_fiber_timeshift.Update();
    }

    cm_fiber_timeshift.ReinitializeHistory(1.2, 0.3);

    for (int timestep = 1; timestep <= 3; ++timestep)
    {
      cm_fiber_timeshift.RecomputeState(1.2 + 0.001 * timestep, 0.3 + 0.1 * timestep);

      ASSERT_DOUBLE_EQ(
          cm_fiber_timeshift.computed_growth_scalar_, cm_fiber_timeshift.computed_growth_scalar_);

      cm_fiber_timeshift.Update();
    }

    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.reference_time_, 0.0);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[0][0].deposition_time, 0.0);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[0][1].deposition_time, 0.1);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[0][2].deposition_time, 0.2);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[0][3].deposition_time, 0.3);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[1][0].deposition_time, 0.3);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[1][1].deposition_time, 0.4);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[1][2].deposition_time, 0.5);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[1][3].deposition_time, 0.6);

    cm_fiber_timeshift.AddTime(1.0);

    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.reference_time_, 1.0);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[0][0].deposition_time, 1.0);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[0][1].deposition_time, 1.1);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[0][2].deposition_time, 1.2);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[0][3].deposition_time, 1.3);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[1][0].deposition_time, 1.3);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[1][1].deposition_time, 1.4);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[1][2].deposition_time, 1.5);
    ASSERT_DOUBLE_EQ(cm_fiber_timeshift.history_[1][3].deposition_time, 1.6);
  }

  TEST_F(FullConstrainedMixtureFiberTest, FiberWithJumpInTimeEqualsWithoutJumpIntime)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber_default = GenerateFiber<double>();
    MIXTURE::FullConstrainedMixtureFiber cm_fiber_timeshift = GenerateFiber<double>();
    std::array all_lambda_f = {1.1, 1.2, 1.3};


    int timestep = 1;
    for (double lambda_f : all_lambda_f)
    {
      cm_fiber_default.ReinitializeHistory(lambda_f, 0.1 * (timestep - 1));
      cm_fiber_timeshift.ReinitializeHistory(lambda_f, 0.1 * (timestep - 1));

      for (int i = 0; i <= 3; ++i)
      {
        cm_fiber_default.RecomputeState(lambda_f + 0.001 * timestep, 0.1 * timestep);
        cm_fiber_timeshift.RecomputeState(lambda_f + 0.001 * timestep, 0.1 * timestep);

        ASSERT_DOUBLE_EQ(cm_fiber_default.computed_sigma_, cm_fiber_default.computed_sigma_);
        ASSERT_DOUBLE_EQ(
            cm_fiber_timeshift.computed_growth_scalar_, cm_fiber_timeshift.computed_growth_scalar_);

        cm_fiber_default.Update();
        cm_fiber_timeshift.Update();

        ++timestep;
      }
    }

    cm_fiber_timeshift.AddTime(1.0);

    cm_fiber_default.RecomputeState(1.3 + 0.001 * timestep, 0.1 * timestep);
    cm_fiber_timeshift.RecomputeState(1.3 + 0.001 * timestep, 0.1 * timestep + 1.0);

    ASSERT_DOUBLE_EQ(cm_fiber_default.computed_sigma_, cm_fiber_default.computed_sigma_);
    ASSERT_DOUBLE_EQ(
        cm_fiber_timeshift.computed_growth_scalar_, cm_fiber_timeshift.computed_growth_scalar_);
  }


  TEST_F(FullConstrainedMixtureFiberTest, GetLastTimeInHistoryIfHistoryIsEmpty)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();

    ASSERT_DOUBLE_EQ(cm_fiber.GetLastTimeInHistory(), 0.0);
  }


  TEST_F(FullConstrainedMixtureFiberTest, GetLastTimeInHistoryIfHistoryIsEmptyAndAfterAddTime)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();

    cm_fiber.AddTime(1.0);
    ASSERT_DOUBLE_EQ(cm_fiber.GetLastTimeInHistory(), 1.0);
  }

  TEST_F(FullConstrainedMixtureFiberTest, GetLastTimeInHistoryWithMultipleIntervals)
  {
    MIXTURE::FullConstrainedMixtureFiber cm_fiber = GenerateFiber<double>();

    cm_fiber.ReinitializeHistory(1.0, 0.0);

    for (int timestep = 1; timestep <= 3; ++timestep)
    {
      cm_fiber.RecomputeState(1.1 + 0.001 * timestep, 0.1 * timestep);

      ASSERT_DOUBLE_EQ(cm_fiber.computed_growth_scalar_, cm_fiber.computed_growth_scalar_);

      cm_fiber.Update();
    }

    cm_fiber.ReinitializeHistory(1.2, 0.3);

    for (int timestep = 1; timestep <= 3; ++timestep)
    {
      cm_fiber.RecomputeState(1.2 + 0.001 * timestep, 0.3 + 0.1 * timestep);

      ASSERT_DOUBLE_EQ(cm_fiber.computed_growth_scalar_, cm_fiber.computed_growth_scalar_);

      cm_fiber.Update();
    }

    ASSERT_DOUBLE_EQ(cm_fiber.GetLastTimeInHistory(), 0.6);
  }

  TEST_F(FullConstrainedMixtureFiberTest, UpdateCallsAddTimeIfGrowthIsDisabled)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>(12.0, 0.1, 1.1, true, false);

    for (int timestep = 1; timestep <= 3; ++timestep)
    {
      const double time = 0.3 + 0.1 * timestep;
      fiber.RecomputeState(1.2 + 0.001 * timestep, time);

      fiber.Update();

      EXPECT_DOUBLE_EQ(fiber.reference_time_, time);
    }
  }

  TEST_F(FullConstrainedMixtureFiberTest,
      CauchyStressAndGrowthScalarIsNormalIfHistoryIsEmptyAndGrowthIsDisabled)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>(12.0, 0.1, 1.1, true, false);
    MIXTURE::FullConstrainedMixtureFiber fiber_nogrowth =
        GenerateFiber<double>(1e100, 0.0, 1.1, false, true);

    fiber_nogrowth.ReinitializeHistory(1.2, 0.0);

    for (int timestep = 1; timestep <= 3; ++timestep)
    {
      fiber.RecomputeState(1.2 + 0.001 * timestep, 0.3 + 0.1 * timestep);
      fiber_nogrowth.RecomputeState(1.2 + 0.001 * timestep, 0.3 + 0.1 * timestep);

      EXPECT_DOUBLE_EQ(fiber.computed_sigma_, fiber_nogrowth.computed_sigma_);
      EXPECT_DOUBLE_EQ(fiber.computed_growth_scalar_, 1.0);
      EXPECT_DOUBLE_EQ(fiber_nogrowth.computed_growth_scalar_, 1.0);

      fiber.Update();
      fiber_nogrowth.Update();
    }
  }

  TEST_F(FullConstrainedMixtureFiberTest,
      CauchyStressAndGrowthScalarRemainNormalIfGrowthIsDisabledAfterGrowthPeriod)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>(12.0, 0.1, 1.1, true, true);

    double lambda_f = 1.2;
    fiber.ReinitializeHistory(lambda_f, 0.0);
    for (int timestep = 1; timestep <= 3; ++timestep)
    {
      lambda_f = 1.2 + 0.001 * timestep;
      fiber.RecomputeState(lambda_f, 0.3 + 0.1 * timestep);

      fiber.Update();
    }

    fiber.growth_enabled_ = false;

    const double growth_scalar = fiber.computed_growth_scalar_;
    const double sigma = fiber.computed_sigma_;

    for (int timestep = 4; timestep <= 6; ++timestep)
    {
      fiber.RecomputeState(lambda_f, 0.3 + 0.1 * timestep);

      EXPECT_DOUBLE_EQ(fiber.computed_growth_scalar_, growth_scalar);
      EXPECT_DOUBLE_EQ(fiber.computed_sigma_, sigma);
      fiber.Update();
    }
  }

  TEST_F(FullConstrainedMixtureFiberTest,
      CauchyStressAndGrowthScalarRemainNormalIfGrowthIsDisabledAfterGrowthPeriodAndInitialGrowthFreePeriod)
  {
    MIXTURE::FullConstrainedMixtureFiber fiber = GenerateFiber<double>(12.0, 0.1, 1.1, true, false);

    double lambda_f = 1.2;

    for (int timestep = 1; timestep <= 3; ++timestep)
    {
      lambda_f = 1.2 + 0.001 * timestep;
      fiber.RecomputeState(lambda_f, 0.3 + 0.1 * timestep);

      EXPECT_DOUBLE_EQ(fiber.computed_growth_scalar_, 1.0);

      fiber.Update();
    }
    fiber.growth_enabled_ = true;

    EXPECT_DOUBLE_EQ(fiber.computed_growth_scalar_, 1.0);

    fiber.ReinitializeHistory(lambda_f, 0.0);
    for (int timestep = 4; timestep <= 6; ++timestep)
    {
      lambda_f = 1.2 + 0.001 * timestep;
      fiber.RecomputeState(lambda_f, 0.3 + 0.1 * timestep);

      fiber.Update();
    }

    fiber.growth_enabled_ = false;

    const double growth_scalar = fiber.computed_growth_scalar_;
    const double sigma = fiber.computed_sigma_;

    for (int timestep = 7; timestep <= 9; ++timestep)
    {
      fiber.RecomputeState(lambda_f, 0.3 + 0.1 * timestep);

      EXPECT_DOUBLE_EQ(fiber.computed_growth_scalar_, growth_scalar);
      EXPECT_DOUBLE_EQ(fiber.computed_sigma_, sigma);
      fiber.Update();
    }
  }

}  // namespace