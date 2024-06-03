/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the ElastHyper service functions
\level 2

*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_elasthyper_service.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_coupanisoexpo.hpp"
#include "4C_matelast_isoneohooke.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace FourC;

  class ElastHyperServiceTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      prinv_(0) = 3.1;
      prinv_(1) = 3.2;
      prinv_(2) = 1.05;
    }

    CORE::LINALG::Matrix<3, 1> prinv_;
  };

  TEST_F(ElastHyperServiceTest, TestCalculateGammaDelta)
  {
    // required inputs to test the method
    // first derivative of strain energy function w.r.t. principle invariants
    CORE::LINALG::Matrix<3, 1> dPI;
    dPI(0) = 0.536;
    dPI(1) = 0.681;
    dPI(2) = 3.981;

    // second derivative of strain energy function w.r.t. principle invariants
    CORE::LINALG::Matrix<6, 1> ddPII;
    ddPII(0) = 0.153;
    ddPII(1) = 0.861;
    ddPII(2) = 2.964;
    ddPII(3) = 1.549;
    ddPII(4) = 0.682;
    ddPII(5) = 3.584;

    // factors for stress calculation as in Holzapfel: Nonlinear Solid Mechanics, 2000:
    // p. 216 (6.32)
    CORE::LINALG::Matrix<3, 1> gamma(true);
    // reference solutions
    CORE::LINALG::Matrix<3, 1> gamma_ref;
    gamma_ref(0) = 5.294200000000001;
    gamma_ref(1) = -1.362;
    gamma_ref(2) = 8.360099999999999;

    // factors for calculation of elasticity tensor as in Holzapfel: Nonlinear Solid Mechanics,
    // 2000: p. 261 (6.194)
    CORE::LINALG::Matrix<8, 1> delta(true);
    // reference solutions
    CORE::LINALG::Matrix<8, 1> delta_ref;
    delta_ref(0) = 125.31604;
    delta_ref(1) = -25.0124;
    delta_ref(2) = 23.03238;
    delta_ref(3) = 3.444;
    delta_ref(4) = -6.5058;
    delta_ref(5) = 29.79144;
    delta_ref(6) = -16.7202;
    delta_ref(7) = -2.724;

    MAT::CalculateGammaDelta(gamma, delta, prinv_, dPI, ddPII);

    FOUR_C_EXPECT_NEAR(gamma, gamma_ref, 1.0e-10);
    FOUR_C_EXPECT_NEAR(delta, delta_ref, 1.0e-10);
  }

  TEST_F(ElastHyperServiceTest, TestEvaluateRightCauchyGreenStrainLikeVoigt)
  {
    // Green-Lagrange tensor in strain-like voigt notation as input
    CORE::LINALG::Matrix<6, 1> E_VoigtStrain;
    E_VoigtStrain(0, 0) = 0.1;
    E_VoigtStrain(1, 0) = 0.2;
    E_VoigtStrain(2, 0) = 0.3;
    E_VoigtStrain(3, 0) = 0.01;
    E_VoigtStrain(4, 0) = 0.02;
    E_VoigtStrain(5, 0) = 0.03;

    CORE::LINALG::Matrix<6, 1> C_VoigtStrain(true);
    // set up reference solution for Cauchy-Green tensor in strain-like voigt notation
    CORE::LINALG::Matrix<6, 1> C_VoigtStrain_ref(true);
    C_VoigtStrain_ref(0) = 1.2;
    C_VoigtStrain_ref(1) = 1.4;
    C_VoigtStrain_ref(2) = 1.6;
    C_VoigtStrain_ref(3) = 0.02;
    C_VoigtStrain_ref(4) = 0.04;
    C_VoigtStrain_ref(5) = 0.06;

    MAT::EvaluateRightCauchyGreenStrainLikeVoigt(E_VoigtStrain, C_VoigtStrain);

    FOUR_C_EXPECT_NEAR(C_VoigtStrain, C_VoigtStrain_ref, 1.0e-10);
  }

  TEST_F(ElastHyperServiceTest, TestEvaluateInvariantDerivatives)
  {
    // Create summands to test the elast hyper functions
    // Currently, only an isotropic summand is supported. Anisotropic summands are not easy to be
    // added here as the setup of the anisotropy is not trivial.

    // Create parameter of IsoNeoHooke material
    IO::InputParameterContainer iso_neo_hooke_data;
    iso_neo_hooke_data.Add("MUE", 1.3);

    auto iso_neo_hooke_params =
        MAT::make_parameter(1, CORE::Materials::MaterialType::mes_isoneohooke, iso_neo_hooke_data);

    // Create summand vector
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsum(0);
    potsum.emplace_back(Teuchos::rcp(new MAT::ELASTIC::IsoNeoHooke(
        dynamic_cast<MAT::ELASTIC::PAR::IsoNeoHooke *>(iso_neo_hooke_params.get()))));

    // Read summand properties
    MAT::SummandProperties properties;
    MAT::ElastHyperProperties(potsum, properties);

    // first derivative of strain energy function w.r.t. principle invariants
    CORE::LINALG::Matrix<3, 1> dPI(true);
    // reference solutions
    CORE::LINALG::Matrix<3, 1> dPI_ref;
    dPI_ref(0) = 0.639514;
    dPI_ref(1) = 0.0;
    dPI_ref(2) = -0.629363;

    // second derivative of strain energy function w.r.t. principle invariants
    CORE::LINALG::Matrix<6, 1> ddPII(true);
    // reference solutions
    CORE::LINALG::Matrix<6, 1> ddPII_ref;
    ddPII_ref(0) = 0.0;
    ddPII_ref(1) = 0.0;
    ddPII_ref(2) = 0.799191;
    ddPII_ref(3) = 0.0;
    ddPII_ref(4) = -0.20302;
    ddPII_ref(5) = 0.0;

    // Compute derivatives of the strain energy function w.r.t. the principal invariants of the
    // strain energy function
    MAT::ElastHyperEvaluateInvariantDerivatives(prinv_, dPI, ddPII, potsum, properties, 0, 0);

    FOUR_C_EXPECT_NEAR(dPI, dPI_ref, 1.0e-4);
    FOUR_C_EXPECT_NEAR(ddPII, ddPII_ref, 1.0e-4);
  }
}  // namespace
