/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the remodel fiber
\level 3
*/
/*----------------------------------------------------------------------*/

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
    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, T> generate_fiber()
    {
      const auto container = Teuchos::rcp(new CORE::MAT::PAR::Material());
      container->Add("K1", 1.3);
      container->Add("K2", 1.3);
      container->Add("COMPRESSION", true);

      fiber_material_parameter_ =
          std::make_shared<MIXTURE::PAR::RemodelFiberMaterialExponential<FADdouble>>(container);

      const auto material =
          std::make_shared<const MIXTURE::RemodelFiberMaterialExponential<FADdouble>>(
              fiber_material_parameter_.get());

      MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, T> fiber(material, {3.4, 12.0}, 1.1);
      return fiber;
    }

    // parameters
    std::shared_ptr<MIXTURE::PAR::RemodelFiberMaterialExponential<FADdouble>>
        fiber_material_parameter_;
  };

  TEST_F(RemodelFiberTest, TestEvaluateDGrowthEvolutionEquationDtDGrowth)
  {
    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, FADdouble> fiber =
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
    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, FADdouble> fiber =
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
    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, FADdouble> fiber =
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
    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, FADdouble> fiber =
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
    MIXTURE::IMPLEMENTATION::RemodelFiberImplementation<2, FADdouble> fiber =
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