// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent_remodelfiber_material_exponential.hpp"
#include "4C_mixture_growth_evolution_linear_cauchy_poisson_turnover.hpp"
#include "4C_mixture_remodelfiber-internal.hpp"
#include "4C_mixture_remodelfiber.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <Sacado.hpp>

#include <memory>

namespace
{
  using namespace FourC;
  using FADdouble = Sacado::Fad::DFad<double>;

  class RemodelFiberTest : public ::testing::Test
  {
   protected:
    template <typename T>
    Mixture::Implementation::RemodelFiberImplementation<2, T> generate_fiber()
    {
      Core::IO::InputParameterContainer container;
      container.add("K1", 1.3);
      container.add("K2", 1.3);
      container.add("COMPRESSION", true);

      fiber_material_parameter_ =
          std::make_shared<Mixture::PAR::RemodelFiberMaterialExponential<FADdouble>>(
              Core::Mat::PAR::Parameter::Data{.parameters = container});

      const auto material =
          std::make_shared<const Mixture::RemodelFiberMaterialExponential<FADdouble>>(
              fiber_material_parameter_.get());

      Mixture::Implementation::RemodelFiberImplementation<2, T> fiber(
          material, {3.4, 12.0, true}, 1.1);
      return fiber;
    }

    // parameters
    std::shared_ptr<Mixture::PAR::RemodelFiberMaterialExponential<FADdouble>>
        fiber_material_parameter_;
  };

  TEST_F(RemodelFiberTest, TestEvaluateDGrowthEvolutionEquationDtDGrowth)
  {
    Mixture::Implementation::RemodelFiberImplementation<2, FADdouble> fiber =
        generate_fiber<FADdouble>();

    const double lambda_f = 1.02;
    const double lambda_ext = 1.014;
    const FADdouble growth_scalar = FADdouble(2, 0, 1.12);
    const FADdouble lambda_r = FADdouble(2, 1, 1.05);

    const FADdouble y =
        fiber.evaluate_growth_evolution_equation_dt(lambda_f, lambda_r, lambda_ext, growth_scalar);
    const FADdouble dGrowthEvolutionEquationDtDGrowth =
        fiber.evaluate_d_growth_evolution_equation_dt_d_growth(
            lambda_f, lambda_r, lambda_ext, growth_scalar);

    EXPECT_FLOAT_EQ(y.dx(0), dGrowthEvolutionEquationDtDGrowth.val());
  }

  TEST_F(RemodelFiberTest, TestEvaluateDGrowthEvolutionEquationDtDRemodel)
  {
    Mixture::Implementation::RemodelFiberImplementation<2, FADdouble> fiber =
        generate_fiber<FADdouble>();

    const double lambda_f = 1.02;
    const double lambda_ext = 1.014;
    const FADdouble growth_scalar = FADdouble(2, 0, 1.12);
    const FADdouble lambda_r = FADdouble(2, 1, 1.05);

    const FADdouble y =
        fiber.evaluate_growth_evolution_equation_dt(lambda_f, lambda_r, lambda_ext, growth_scalar);
    const FADdouble dGrowthEvolutionEquationDtDRemodel =
        fiber.evaluate_d_growth_evolution_equation_dt_d_remodel(
            lambda_f, lambda_r, lambda_ext, growth_scalar);

    EXPECT_FLOAT_EQ(y.dx(1), dGrowthEvolutionEquationDtDRemodel.val());
  }

  TEST_F(RemodelFiberTest, TestEvaluateDRemodelEvolutionEquationDtDGrowth)
  {
    Mixture::Implementation::RemodelFiberImplementation<2, FADdouble> fiber =
        generate_fiber<FADdouble>();

    const double lambda_f = 1.02;
    const double lambda_ext = 1.014;
    const FADdouble growth_scalar = FADdouble(2, 0, 1.12);
    const FADdouble lambda_r = FADdouble(2, 1, 1.05);

    const FADdouble y =
        fiber.evaluate_remodel_evolution_equation_dt(lambda_f, lambda_r, lambda_ext);
    const FADdouble dRemodelEvolutionEquationDtDGrowth =
        fiber.evaluate_d_remodel_evolution_equation_dt_d_growth(lambda_f, lambda_r, lambda_ext);

    EXPECT_FLOAT_EQ(y.dx(0), dRemodelEvolutionEquationDtDGrowth.val());
  }

  TEST_F(RemodelFiberTest, TestEvaluateDRemodelEvolutionEquationDtDRemodel)
  {
    Mixture::Implementation::RemodelFiberImplementation<2, FADdouble> fiber =
        generate_fiber<FADdouble>();

    const double lambda_f = 1.02;
    const double lambda_ext = 1.014;
    const FADdouble growth_scalar = FADdouble(2, 0, 1.12);
    const FADdouble lambda_r = FADdouble(2, 1, 1.05);

    const FADdouble y =
        fiber.evaluate_remodel_evolution_equation_dt(lambda_f, lambda_r, lambda_ext);
    const FADdouble dRemodelEvolutionEquationDtDRemodel =
        fiber.evaluate_d_remodel_evolution_equation_dt_d_remodel(lambda_f, lambda_r, lambda_ext);

    EXPECT_FLOAT_EQ(y.dx(1), dRemodelEvolutionEquationDtDRemodel.val());
  }

  TEST_F(RemodelFiberTest, TestEvaluateDFiberCauchyStressDRemodel)
  {
    Mixture::Implementation::RemodelFiberImplementation<2, FADdouble> fiber =
        generate_fiber<FADdouble>();

    const double lambda_f = 1.02;
    const double lambda_ext = 1.014;
    const FADdouble growth_scalar = FADdouble(2, 0, 1.12);
    const FADdouble lambda_r = FADdouble(2, 1, 1.05);

    const FADdouble y = fiber.evaluate_fiber_cauchy_stress(lambda_f, lambda_r, lambda_ext);
    const FADdouble dFiberCauchyStressDRemodel =
        fiber.evaluate_d_fiber_cauchy_stress_d_remodel(lambda_f, lambda_r, lambda_ext);

    EXPECT_FLOAT_EQ(y.dx(1), dFiberCauchyStressDRemodel.val());
  }
}  // namespace
