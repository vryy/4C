/*----------------------------------------------------------------------*/
/*! \file
\brief unit testing functionality for druckerprager material
\level 3
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_beam3_kirchhoff.hpp"
#include "4C_beam3_reissner.hpp"
#include "4C_comm_pack_buffer.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_FADmatrix_utils.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_plasticdruckerprager.hpp"
#include "4C_mat_service.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
namespace
{
  using namespace FourC;

  class DruckerPragerTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      Core::IO::InputParameterContainer container;
      container.Add("YOUNG", 1.0);
      container.Add("NUE", 0.25);
      container.Add("DENS", 0.0);
      container.Add("ISOHARD", 1.0);
      container.Add("TOL", 1.e-12);
      container.Add("C", 1.);
      container.Add("ETA", 1.);
      container.Add("XI", 1.);
      container.Add("ETABAR", 1.);
      container.Add("MAXITER", 50);
      container.Add("TANG", std::string("consistent"));

      param_druckprag_ = std::shared_ptr(
          Mat::make_parameter(1, Core::Materials::MaterialType::m_pldruckprag, container));

      Global::Problem& problem = (*Global::Problem::Instance());
      problem.Materials()->SetReadFromProblem(0);
      problem.Materials()->insert(1, param_druckprag_);
      problem.Materials().assert_not_null();
      druckprag_ = Teuchos::rcp(new Mat::PlasticDruckerPrager(
          dynamic_cast<Mat::PAR::PlasticDruckerPrager*>(param_druckprag_.get())));
    }

    void TearDown() override
    {
      // We need to make sure the Global::Problem instance created in SetUp is deleted again. If
      // this is not done, some troubles arise where unit tests influence each other on some
      // configurations. We suspect that missing singleton destruction might be the reason for that.
      Global::Problem::Done();
    };
    std::shared_ptr<Core::Mat::PAR::Parameter> param_druckprag_;
    Core::Communication::PackBuffer data;
    Teuchos::RCP<Mat::PlasticDruckerPrager> druckprag_;
  };

  //! test member function Pack and unpack
  TEST_F(DruckerPragerTest, TestPackUnpack)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1> input_strain;
    for (int i = 0; i < 6; ++i) input_strain(i) = .1;
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    Core::LinAlg::Matrix<6, 1> ref_stress(true);
    for (int i = 0; i < 3; ++i)
      ref_stress(i) =
          (1.0 / ((1.0 + 0.25) * (1.0 - (2.0 * 0.25)))) * ((1.0 - 0.25) + 0.25 + 0.25) * .1;
    for (int i = 3; i < 6; ++i)
      ref_stress(i) =
          (1.0 / ((1.0 + 0.25) * (1.0 - (2.0 * 0.25)))) * ((1.0 - (2.0 * 0.25)) / 2.0) * .1;
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1> result_stress(true);
    druckprag_->pack(data);
    std::vector<char> dataSend;
    swap(dataSend, data());
    for (int i = 0; i < 4; i++) dataSend.erase(dataSend.begin());
    auto plastic = Teuchos::rcp(new Mat::PlasticDruckerPrager());
    plastic->unpack(dataSend);
    plastic->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    FOUR_C_EXPECT_NEAR(result_stress, ref_stress, 1.0e-12);
  };

  //! test member function Evaluate
  TEST_F(DruckerPragerTest, TestEvaluate)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1> input_strain;
    for (int i = 0; i < 6; ++i) input_strain(i) = .1;
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    Core::LinAlg::Matrix<6, 1> ref_stress(true);
    for (int i = 0; i < 3; ++i)
      ref_stress(i) =
          (1.0 / ((1.0 + 0.25) * (1.0 - (2.0 * 0.25)))) * ((1.0 - 0.25) + 0.25 + 0.25) * .1;
    for (int i = 3; i < 6; ++i)
      ref_stress(i) =
          (1.0 / ((1.0 + 0.25) * (1.0 - (2.0 * 0.25)))) * ((1.0 - (2.0 * 0.25)) / 2.0) * .1;
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1> result_stress(true);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    FOUR_C_EXPECT_NEAR(result_stress, ref_stress, 1.0e-12);
  };

  //! test member function Evaluate for Return to Cone
  TEST_F(DruckerPragerTest, TestEvaluateReturnToCone)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1> input_strain;
    for (int i = 0; i < 3; ++i) input_strain(i) = 0.0;
    for (int i = 3; i < 6; ++i) input_strain(i) = 2.2;
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    double Dgamma = (2.2 * sqrt(3) / 2.5 - 1.0) / 31.0 * 15.0;
    Core::LinAlg::Matrix<6, 1> ref_stress;
    for (int i = 0; i < 3; ++i) ref_stress(i) = -(Dgamma * (1.0 / (3.0 * (1.0 - (2.0 * 0.25)))));
    for (int i = 3; i < 6; ++i)
      ref_stress(i) = (1.0 / (2 * (1.0 + 0.25))) *
                      (1 - ((1.0 / (2 * (1.0 + 0.25))) * Dgamma / (2.2 * sqrt(3) / 2.5))) * 2.2;
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1> result_stress(true);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    FOUR_C_EXPECT_NEAR(result_stress, ref_stress, 1.0e-12);
  };

  //! test member function Evaluate for Return to Apex
  TEST_F(DruckerPragerTest, TestEvaluateReturnToApex)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1> input_strain;
    for (int i = 0; i < 3; ++i) input_strain(i) = 1.0;
    for (int i = 3; i < 6; ++i) input_strain(i) = 0.0;
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    Core::LinAlg::Matrix<6, 1> ref_stress(true);
    for (int i = 0; i < 3; ++i) ref_stress(i) = 2.0 - (10. / 15.) * (3. / 5.);
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1> result_stress(true);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    FOUR_C_EXPECT_NEAR(result_stress, ref_stress, 1.0e-12);
  };

  //! test member function Evaluate for History and elastic unloading
  TEST_F(DruckerPragerTest, TestEvaluateHistory)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1, FAD> input_strain;
    for (int i = 0; i < 3; ++i) input_strain(i) = FAD(6, i, 0.1);
    for (int i = 3; i < 6; ++i) input_strain(i) = FAD(6, i, 0.1);
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    Core::LinAlg::Matrix<6, 1, FAD> ref_stress(true);
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1, FAD> result_stress(true);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    Core::LinAlg::Matrix<6, 6> ref_cmat(true);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        ref_cmat(i, j) = result_stress(i).dx(j);
      }
    }
    FOUR_C_EXPECT_NEAR(result_cmat, ref_cmat, 1.0e-12);
    druckprag_->update();
    for (int i = 0; i < 3; ++i) input_strain(i) = FAD(6, i, 1.0);
    for (int i = 3; i < 6; ++i) input_strain(i) = FAD(6, i, 0.0);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        ref_cmat(i, j) = result_stress(i).dx(j);
      }
    }
    FOUR_C_EXPECT_NEAR(result_cmat, ref_cmat, 1.0e-12);
    druckprag_->update();
    for (int i = 0; i < 3; ++i) input_strain(i) = FAD(6, i, 0.2);
    for (int i = 3; i < 6; ++i) input_strain(i) = FAD(6, i, 0.0);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        ref_cmat(i, j) = result_stress(i).dx(j);
      }
    }
    FOUR_C_EXPECT_NEAR(Core::FADUtils::CastToDouble(result_stress),
        Core::FADUtils::CastToDouble(ref_stress), 1.0e-12);
  };

  //! test member function Evaluate for arbitrary values
  TEST_F(DruckerPragerTest, TestEvaluateRandomStrain)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1> input_strain;
    input_strain(0) = 1.1;
    input_strain(1) = 2.0;
    input_strain(2) = 0.1;
    input_strain(3) = 2.5;
    input_strain(4) = 1.4;
    input_strain(5) = 1.0;
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    Core::LinAlg::Matrix<6, 1> ref_stress(true);
    ref_stress(0) = 1.3231031817668;
    ref_stress(1) = 1.7934880206154;
    ref_stress(2) = 0.8004533608238;
    ref_stress(3) = 0.6533122761787;
    ref_stress(4) = 0.3658548746601;
    ref_stress(5) = 0.2613249104715;
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1> result_stress(true);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    FOUR_C_EXPECT_NEAR(result_stress, ref_stress, 1.0e-12);
  };

  //! test member function Evaluate
  TEST_F(DruckerPragerTest, TestEvaluateCmat)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1, FAD> input_strain;
    for (int i = 0; i < 6; ++i) input_strain(i) = FAD(6, i, .1 * i);
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    Core::LinAlg::Matrix<6, 1, FAD> ref_stress(true);
    for (int i = 0; i < 3; ++i)
      ref_stress(i) =
          FAD((1.0 / ((1.0 + 0.25) * (1.0 - (2.0 * 0.25)))) * ((1.0 - 0.25) + 0.25 + 0.25) * .1);
    for (int i = 3; i < 6; ++i)
      ref_stress(i) =
          FAD((1.0 / ((1.0 + 0.25) * (1.0 - (2.0 * 0.25)))) * ((1.0 - (2.0 * 0.25)) / 2.0) * .1);
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1, FAD> result_stress(true);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    Core::LinAlg::Matrix<6, 6> ref_cmat(true);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        ref_cmat(i, j) = result_stress(i).dx(j);
      }
    }
    FOUR_C_EXPECT_NEAR(result_cmat, ref_cmat, 1.0e-12);
  };

  //! test CMAT matrix for Return to Cone
  TEST_F(DruckerPragerTest, TestEvaluateReturnToConeCmat)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1, FAD> input_strain;
    for (int i = 0; i < 3; ++i) input_strain(i) = FAD(6, i, 0.1 * i);
    for (int i = 3; i < 6; ++i) input_strain(i) = FAD(6, i, 2.2 * i);
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1, FAD> result_stress(true);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    Core::LinAlg::Matrix<6, 6> ref_cmat(true);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        ref_cmat(i, j) = result_stress(i).dx(j);
      }
    }
    FOUR_C_EXPECT_NEAR(result_cmat, ref_cmat, 1.0e-12);
  };
  TEST_F(DruckerPragerTest, TestEvaluateReturnToApexCmat)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1, FAD> input_strain;
    for (int i = 0; i < 3; ++i) input_strain(i) = FAD(6, i, 1.0);
    for (int i = 3; i < 6; ++i) input_strain(i) = FAD(6, i, 0.0);
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1, FAD> result_stress(true);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    Core::LinAlg::Matrix<6, 6> ref_cmat(true);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        ref_cmat(i, j) = result_stress(i).dx(j);
      }
    }
    FOUR_C_EXPECT_NEAR(result_cmat, ref_cmat, 1.0e-12);
  };

  //! test CMAT matrix for Return to Apex
  TEST_F(DruckerPragerTest, TestEvaluateRandomStrainCmat)
  {
    Input::LineDefinition linedef;
    druckprag_->setup(1, &linedef);
    Core::LinAlg::Matrix<6, 1, FAD> input_strain;
    input_strain(0) = FAD(6, 0, 1.1);
    input_strain(1) = FAD(6, 1, 2.0);
    input_strain(2) = FAD(6, 2, 0.1);
    input_strain(3) = FAD(6, 3, 2.5);
    input_strain(4) = FAD(6, 4, 1.4);
    input_strain(5) = FAD(6, 5, 1.0);
    Teuchos::ParameterList paras;
    Core::LinAlg::Matrix<3, 3> defgrad(true);
    Core::LinAlg::Matrix<6, 1, FAD> ref_stress(true);
    ref_stress(0) = FAD(1.4142412329012);
    ref_stress(1) = FAD(1.8571566160540);
    ref_stress(2) = FAD(0.9221130293981);
    ref_stress(3) = FAD(0.6151602543789);
    ref_stress(4) = FAD(0.3444897424522);
    ref_stress(5) = FAD(0.2460641017516);
    Core::LinAlg::Matrix<6, 6> result_cmat(true);
    Core::LinAlg::Matrix<6, 1, FAD> result_stress(true);
    druckprag_->evaluate(&defgrad, &input_strain, paras, &result_stress, &result_cmat, 0, 0);
    Core::LinAlg::Matrix<6, 6> ref_cmat(true);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        ref_cmat(i, j) = result_stress(i).dx(j);
      }
    }
    FOUR_C_EXPECT_NEAR(result_cmat, ref_cmat, 1.0e-12);
  };
}  // namespace
