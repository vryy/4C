
/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the CoupAnisoExpoBase summand with gauss point fibers

\level 2


*/
/*----------------------------------------------------------------------*/


#include <gtest/gtest.h>

#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_fixedsizematrix_voigt_notation.hpp"
#include "baci_mat_anisotropy.hpp"
#include "baci_matelast_aniso_structuraltensor_strategy.hpp"
#include "baci_matelast_coupanisoexpo.hpp"
#include "baci_unittest_utils_assertions_test.hpp"

#include <boost/mpl/protect.hpp>
#include <Teuchos_RCPDecl.hpp>



namespace
{
  using namespace BACI;

  class CoupAnisoExpoAnisotropyExtensionGaussPointFiberTest
      : public ::testing::TestWithParam<std::tuple<int, int>>
  {
   protected:
    CoupAnisoExpoAnisotropyExtensionGaussPointFiberTest()
        : anisotropy_(),
          gpFibers_(2, std::vector<CORE::LINALG::Matrix<3, 1>>(2)),
          gpTensors_(2, std::vector<CORE::LINALG::Matrix<3, 3>>(2)),
          gpTensors_stress_(2, std::vector<CORE::LINALG::Matrix<6, 1>>(2))
    {
      /// initialize dummy fibers
      // gp 0
      gpFibers_[0][0](0) = 0.469809238649817;
      gpFibers_[0][0](1) = 0.872502871778232;
      gpFibers_[0][0](2) = 0.134231211042805;

      gpFibers_[0][1](0) = 0.071428571428571;
      gpFibers_[0][1](1) = 0.142857142857143;
      gpFibers_[0][1](2) = 0.214285714285714;

      // gp 1
      gpFibers_[1][0](0) = 0.245358246032859;
      gpFibers_[1][0](1) = 0.858753861115007;
      gpFibers_[1][0](2) = 0.449823451060242;

      gpFibers_[1][1](0) = 0.068965517241379;
      gpFibers_[1][1](1) = 0.103448275862069;
      gpFibers_[1][1](2) = 0.137931034482759;

      for (std::size_t gp = 0; gp < 2; ++gp)
      {
        for (std::size_t i = 0; i < 2; ++i)
        {
          gpTensors_[gp][i].MultiplyNT(gpFibers_[gp][i], gpFibers_[gp][i]);
          CORE::LINALG::VOIGT::Stresses::MatrixToVector(
              gpTensors_[gp][i], gpTensors_stress_[gp][i]);
        }
      }

      SetupAnisotropyExtension();
    }

    void SetupAnisotropyExtension()
    {
      int fiber_id = std::get<0>(GetParam());
      auto strategy = Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr));
      anisotropyExtension_ = std::make_unique<MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension>(
          3, 0.0, false, strategy, fiber_id);
      anisotropyExtension_->RegisterNeededTensors(
          MAT::FiberAnisotropyExtension<1>::FIBER_VECTORS |
          MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
          MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
      anisotropy_.RegisterAnisotropyExtension(*anisotropyExtension_);
      anisotropy_.SetNumberOfGaussPoints(2);

      // Setup Gauss point fibers
      anisotropy_.SetGaussPointFibers(gpFibers_);
    }

    [[nodiscard]] int GetGaussPoint() const { return std::get<1>(GetParam()); }

    [[nodiscard]] int GetFiberId() const { return std::get<0>(GetParam()); }

    MAT::Anisotropy anisotropy_;
    std::unique_ptr<MAT::ELASTIC::CoupAnisoExpoAnisotropyExtension> anisotropyExtension_;

    std::vector<std::vector<CORE::LINALG::Matrix<3, 1>>> gpFibers_;
    std::vector<std::vector<CORE::LINALG::Matrix<3, 3>>> gpTensors_;
    std::vector<std::vector<CORE::LINALG::Matrix<6, 1>>> gpTensors_stress_;
  };

  TEST_P(CoupAnisoExpoAnisotropyExtensionGaussPointFiberTest, GetScalarProduct)
  {
    EXPECT_NEAR(anisotropyExtension_->GetScalarProduct(GetGaussPoint()), 1.0, 1e-10);
  }

  TEST_P(CoupAnisoExpoAnisotropyExtensionGaussPointFiberTest, GetFiber)
  {
    BACI_EXPECT_NEAR(anisotropyExtension_->GetFiber(GetGaussPoint()),
        gpFibers_[GetGaussPoint()][GetFiberId() - 1], 1e-10);
  }

  TEST_P(CoupAnisoExpoAnisotropyExtensionGaussPointFiberTest, GetStructuralTensor)
  {
    BACI_EXPECT_NEAR(anisotropyExtension_->GetStructuralTensor(GetGaussPoint()),
        gpTensors_[GetGaussPoint()][GetFiberId() - 1], 1e-10);
  }

  TEST_P(CoupAnisoExpoAnisotropyExtensionGaussPointFiberTest, GetStructuralTensorStress)
  {
    BACI_EXPECT_NEAR(anisotropyExtension_->GetStructuralTensor_stress(GetGaussPoint()),
        gpTensors_stress_[GetGaussPoint()][GetFiberId() - 1], 1e-10);
  }

  INSTANTIATE_TEST_SUITE_P(GaussPoints, CoupAnisoExpoAnisotropyExtensionGaussPointFiberTest,
      ::testing::Combine(::testing::Values(1, 2), ::testing::Values(0, 1)));
}  // namespace