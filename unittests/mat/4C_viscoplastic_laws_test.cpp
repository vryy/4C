// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_viscoplastic_laws.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <memory>
namespace
{
  using namespace FourC;

  class ViscoplasticLawsTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      // set up equivalent stress and plastic strain
      equiv_stress_ = 1000.0;
      equiv_plastic_strain_ = 0.001;

      // manually create viscoplastic law (ViscoplasticLawReformulatedJohnsonCook)
      Core::IO::InputParameterContainer viscoplastic_law_reformulated_JC_data;
      viscoplastic_law_reformulated_JC_data.add("STRAIN_RATE_PREFAC", 1.0);
      viscoplastic_law_reformulated_JC_data.add("STRAIN_RATE_EXP_FAC", 0.014);
      viscoplastic_law_reformulated_JC_data.add("INIT_YIELD_STRENGTH", 792.0);
      viscoplastic_law_reformulated_JC_data.add("ISOTROP_HARDEN_PREFAC", 510.0);
      viscoplastic_law_reformulated_JC_data.add("ISOTROP_HARDEN_EXP", 0.26);
      params_viscoplastic_law_reformulated_JC_ =
          std::dynamic_pointer_cast<Mat::PAR::ViscoplasticLawReformulatedJohnsonCook>(
              std::shared_ptr(Mat::make_parameter(1,
                  Core::Materials::MaterialType::mvl_reformulated_Johnson_Cook,
                  viscoplastic_law_reformulated_JC_data)));
      viscoplastic_law_reformulated_JC_ =
          std::make_shared<Mat::ViscoplasticLawReformulatedJohnsonCook>(
              params_viscoplastic_law_reformulated_JC_.get());

      // define setup parameter for InelasticDefGradTransvIsotropElastViscoplast
      Core::IO::InputParameterContainer setup_viscoplastic_law_reformulated_JC;  // can stay empty

      // call setup method for ViscoplasticLawReformulatedJohnsonCook
      int numgp = 8;  // HEX8 element, although not really relevant for the tested methods
      viscoplastic_law_reformulated_JC_->setup(numgp, setup_viscoplastic_law_reformulated_JC);

      // call pre_evaluate
      viscoplastic_law_reformulated_JC_->pre_evaluate(0);
    }

    // equivalent stress
    double equiv_stress_;
    // equivalent stress
    double equiv_plastic_strain_;
    // reference solution for the stress ratio (ViscoplasticLawReformulatedJohnsonCook)
    double stress_ratio_reformulated_JC_solution_;
    // reference solution for the plastic strain rate (ViscoplasticLawReformulatedJohnsonCook)
    double plastic_strain_rate_reformulated_JC_solution_;
    // reference solution for the plastic strain rate derivatives, w.r.t. equivalent stress and
    // plastic strain (ViscoplasticLawReformulatedJohnsonCook)
    Core::LinAlg::Matrix<2, 1> deriv_plastic_strain_rate_reformulated_JC_solution_;
    // pointer to ViscoplasticLawReformulatedJohnsonCook
    std::shared_ptr<Mat::ViscoplasticLawReformulatedJohnsonCook> viscoplastic_law_reformulated_JC_;
    // pointer to parameters of ViscoplasticLawReformulatedJohnsonCook
    std::shared_ptr<Mat::PAR::ViscoplasticLawReformulatedJohnsonCook>
        params_viscoplastic_law_reformulated_JC_;

    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard;
  };

  TEST_F(ViscoplasticLawsTest, TestEvaluateStressRatio)
  {
    // set reference solution
    stress_ratio_reformulated_JC_solution_ = 1.14072049868917;

    // compute solution from the viscoplasticity law
    double stress_ratio_reformulated_JC = viscoplastic_law_reformulated_JC_->evaluate_stress_ratio(
        equiv_stress_, equiv_plastic_strain_);

    // compare solutions
    EXPECT_NEAR(stress_ratio_reformulated_JC_solution_, stress_ratio_reformulated_JC, 1.0e-8);
  }

  TEST_F(ViscoplasticLawsTest, TestEvaluatePlasticStrainRate)
  {
    // set reference solution
    plastic_strain_rate_reformulated_JC_solution_ = 23188.7161986626;

    // compute solution from the viscoplasticity law
    double plastic_strain_rate_reformulated_JC =
        viscoplastic_law_reformulated_JC_->evaluate_plastic_strain_rate(
            equiv_stress_, equiv_plastic_strain_, 0.0, false);

    // compare solutions
    EXPECT_NEAR(
        plastic_strain_rate_reformulated_JC_solution_, plastic_strain_rate_reformulated_JC, 1.0e-8);
  }

  TEST_F(ViscoplasticLawsTest, TestEvaluatePlasticStrainRateDerivatives)
  {
    // set reference solution
    deriv_plastic_strain_rate_reformulated_JC_solution_(0, 0) = 1889.49890189991;
    deriv_plastic_strain_rate_reformulated_JC_solution_(1, 0) = -47431778.9968811;

    // compute solution from the viscoplasticity law
    Core::LinAlg::Matrix<2, 1> deriv_plastic_strain_rate_reformulated_JC =
        viscoplastic_law_reformulated_JC_->evaluate_derivatives_of_plastic_strain_rate(
            equiv_stress_, equiv_plastic_strain_, 0.0, false);

    // compare solutions
    FOUR_C_EXPECT_NEAR(deriv_plastic_strain_rate_reformulated_JC_solution_,
        deriv_plastic_strain_rate_reformulated_JC, 1.0e-6);
  }

}  // namespace
