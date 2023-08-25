/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the multiplicative_split_defgrad_elasthyper_service functions
\level 3

*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_linalg_fixedsizematrix.H"
#include "baci_mat_multiplicative_split_defgrad_elasthyper_service.H"
#include "baci_mat_par_material.H"
#include "baci_matelast_isoneohooke.H"
#include "baci_unittest_utils_assertions_test.H"

namespace
{
  class MultiplicativeSplitDefgradElastHyperServiceTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      FM_(0, 0) = 1.1;
      FM_(1, 1) = 1.2;
      FM_(2, 2) = 1.3;
      FM_(0, 1) = FM_(1, 0) = 0.01;
      FM_(1, 2) = FM_(2, 1) = 0.02;
      FM_(0, 2) = FM_(2, 0) = 0.03;

      iFinM_(0, 0) = 1.04;
      iFinM_(1, 1) = 1.03;
      iFinM_(2, 2) = 1.02;
      iFinM_(0, 1) = iFinM_(1, 0) = 0.003;
      iFinM_(1, 2) = iFinM_(2, 1) = 0.001;
      iFinM_(0, 2) = iFinM_(2, 0) = 0.005;

      CM_.MultiplyTN(FM_, FM_);
      iCinM_.MultiplyNT(iFinM_, iFinM_);
    }

    CORE::LINALG::Matrix<3, 3> FM_;
    CORE::LINALG::Matrix<3, 3> iFinM_;

    CORE::LINALG::Matrix<3, 3> CM_;
    CORE::LINALG::Matrix<3, 3> iCinM_;
  };

  TEST_F(MultiplicativeSplitDefgradElastHyperServiceTest, TestEvaluateCe)
  {
    CORE::LINALG::Matrix<3, 3> CeM_target(false);
    CeM_target(0, 0) = 1.3107725000000006;
    CeM_target(1, 1) = 1.5284889394999996;
    CeM_target(2, 2) = 1.7604995235000003;
    CeM_target(0, 1) = CeM_target(1, 0) = 0.03385382080000001;
    CeM_target(0, 2) = CeM_target(2, 0) = 0.091697784;
    CeM_target(1, 2) = CeM_target(2, 1) = 0.0564151401;

    CORE::LINALG::Matrix<3, 3> CeM(false);
    MAT::EvaluateCe(FM_, iFinM_, CeM);

    BACI_EXPECT_NEAR(CeM, CeM_target, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperServiceTest, TestEvaluateiCinCiCin)
  {
    CORE::LINALG::Matrix<3, 3> iCinCiCinM_target(false);
    iCinCiCinM_target(0, 0) = 1.418955902138138;
    iCinCiCinM_target(1, 1) = 1.6219134553554275;
    iCinCiCinM_target(2, 2) = 1.832708744871652;
    iCinCiCinM_target(0, 1) = iCinCiCinM_target(1, 0) = 0.045473409425074995;
    iCinCiCinM_target(0, 2) = iCinCiCinM_target(2, 0) = 0.113283079933819;
    iCinCiCinM_target(1, 2) = iCinCiCinM_target(2, 1) = 0.0631150197598975;

    CORE::LINALG::Matrix<3, 3> iCinCiCinM(false);
    MAT::EvaluateiCinCiCin(CM_, iCinM_, iCinCiCinM);

    BACI_EXPECT_NEAR(iCinCiCinM, iCinCiCinM_target, 1.0e-10);
  }


  TEST_F(MultiplicativeSplitDefgradElastHyperServiceTest, TestElastHyperEvaluateElasticPart)
  {
    CORE::LINALG::Matrix<6, 1> S_stress;
    CORE::LINALG::Matrix<6, 6> cmat;

    // Create parameter of IsoNeoHooke material
    Teuchos::RCP<MAT::PAR::Material> isoNeoHookeParams = Teuchos::rcp(new MAT::PAR::Material());
    isoNeoHookeParams->Add("MUE", 1.3);
    isoNeoHookeParams->SetParameter(new MAT::ELASTIC::PAR::IsoNeoHooke(isoNeoHookeParams));
    auto* isoNeoHookeParams2 =
        dynamic_cast<MAT::ELASTIC::PAR::IsoNeoHooke*>(isoNeoHookeParams->Parameter());

    // Create summand vector
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsum(0);
    potsum.emplace_back(Teuchos::rcp(new MAT::ELASTIC::IsoNeoHooke(isoNeoHookeParams2)));

    // Read summand properties
    MAT::SummandProperties properties;
    MAT::ElastHyperProperties(potsum, properties);

    // Evaluate method to test
    MAT::ElastHyperEvaluateElasticPart(FM_, iFinM_, S_stress, cmat, potsum, properties, 0, 0);

    // Build matrices with the correct solution
    CORE::LINALG::Matrix<6, 1> S_stress_target;
    CORE::LINALG::Matrix<6, 6> cmat_target;

    S_stress_target(0) = -0.16087311035295149;
    S_stress_target(1) = -0.0041652481745947378;
    S_stress_target(2) = 0.11178774018437565;
    S_stress_target(3) = 0.021510980954630672;
    S_stress_target(4) = 0.028193006189638534;
    S_stress_target(5) = 0.054703856893060218;
    cmat_target(0, 0) = 1.3769715241022054;
    cmat_target(0, 1) = -0.42612503282604719;
    cmat_target(0, 2) = -0.42455257184204381;
    cmat_target(0, 3) = -0.031079302592187129;
    cmat_target(0, 4) = -0.0020029032555127484;
    cmat_target(0, 5) = -0.084619011701144234;
    cmat_target(1, 0) = -0.42612503282604719;
    cmat_target(1, 1) = 0.84893700719517029;
    cmat_target(1, 2) = -0.40977461472058613;
    cmat_target(1, 3) = -0.025011243966307363;
    cmat_target(1, 4) = -0.037612994171778887;
    cmat_target(1, 5) = -0.0035293902859759667;
    cmat_target(2, 0) = -0.42455257184204381;
    cmat_target(2, 1) = -0.40977461472058613;
    cmat_target(2, 2) = 0.52732888530268807;
    cmat_target(2, 3) = 0.0006656868830938736;
    cmat_target(2, 4) = -0.030557790648953649;
    cmat_target(2, 5) = -0.055336798871663345;
    cmat_target(3, 0) = -0.031079302592187129;
    cmat_target(3, 1) = -0.025011243966307363;
    cmat_target(3, 2) = 0.0006656868830938736;
    cmat_target(3, 3) = 0.75548157994492326;
    cmat_target(3, 4) = -0.031268379573823023;
    cmat_target(3, 5) = -0.020742640930106387;
    cmat_target(4, 0) = -0.0020029032555127484;
    cmat_target(4, 1) = -0.037612994171778887;
    cmat_target(4, 2) = -0.030557790648953649;
    cmat_target(4, 3) = -0.031268379573823023;
    cmat_target(4, 4) = 0.54196060823048886;
    cmat_target(4, 5) = -0.0079082879496605533;
    cmat_target(5, 0) = -0.084619011701144234;
    cmat_target(5, 1) = -0.0035293902859759667;
    cmat_target(5, 2) = -0.055336798871663345;
    cmat_target(5, 3) = -0.020742640930106387;
    cmat_target(5, 4) = -0.0079082879496605533;
    cmat_target(5, 5) = 0.64761600029220823;

    BACI_EXPECT_NEAR(S_stress, S_stress_target, 1.0e-9);
    BACI_EXPECT_NEAR(cmat, cmat_target, 1.0e-9);
  }
}  // namespace
