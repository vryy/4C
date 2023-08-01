/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the CoupAnisoExpoShear summand with gauss point fibers

\level 2


*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_mat_anisotropy.H"
#include "baci_matelast_coupanisoexpo.H"
#include "baci_lib_voigt_notation.H"
#include "baci_matelast_coupanisoexposhear.H"
#include "baci_unittest_utils_assertions.h"

namespace
{


  void SetupSingleStructuralTensor(const CORE::LINALG::Matrix<3, 1>& fiber1,
      const CORE::LINALG::Matrix<3, 1>& fiber2, CORE::LINALG::Matrix<3, 3>& structuralTensor)
  {
    CORE::LINALG::Matrix<3, 3> fiber1fiber2T(false);

    fiber1fiber2T.MultiplyNT(fiber1, fiber2);

    structuralTensor.Update(0.5, fiber1fiber2T);
    structuralTensor.UpdateT(0.5, fiber1fiber2T, 1.0);
  }

  class CoupAnisoExpoShearGaussPointFibersTest
      : public ::testing::TestWithParam<std::tuple<std::array<int, 2>, int>>
  {
   public:
    CoupAnisoExpoShearGaussPointFibersTest()
        : anisotropy_(),
          gpFibers_(2, std::vector<CORE::LINALG::Matrix<3, 1>>(3)),
          gpTensors_(2),
          gpTensors_stress_(2),
          gpScalarProducts_(2)
    {
      /// initialize dummy fibers
      // gp 0
      gpFibers_[0][0](0) = 0.469809238649817;
      gpFibers_[0][0](1) = 0.872502871778232;
      gpFibers_[0][0](2) = 0.134231211042805;

      gpFibers_[0][1](0) = 0.071428571428571;
      gpFibers_[0][1](1) = 0.142857142857143;
      gpFibers_[0][1](2) = 0.214285714285714;

      gpFibers_[0][2](0) = 0.142857142857143;
      gpFibers_[0][2](1) = 0.214285714285714;
      gpFibers_[0][2](2) = 0.071428571428571;

      // gp 1
      gpFibers_[1][0](0) = 0.245358246032859;
      gpFibers_[1][0](1) = 0.858753861115007;
      gpFibers_[1][0](2) = 0.449823451060242;

      gpFibers_[1][1](0) = 0.068965517241379;
      gpFibers_[1][1](1) = 0.103448275862069;
      gpFibers_[1][1](2) = 0.137931034482759;

      gpFibers_[1][2](0) = 0.137931034482759;
      gpFibers_[1][2](1) = 0.068965517241379;
      gpFibers_[1][2](2) = 0.103448275862069;

      for (std::size_t gp = 0; gp < 2; ++gp)
      {
        // setup structural tensor
        SetupSingleStructuralTensor(
            gpFibers_[gp][GetFiberIds()[0]], gpFibers_[gp][GetFiberIds()[1]], gpTensors_[gp]);

        // Setup structural tensors in stress like Voigt notation
        UTILS::VOIGT::Stresses::MatrixToVector(gpTensors_[gp], gpTensors_stress_[gp]);

        // setup scalar product
        gpScalarProducts_[gp] =
            gpFibers_[gp][GetFiberIds()[0]].Dot(gpFibers_[gp][GetFiberIds()[1]]);
      }

      SetupAnisotropyExtension(GetFiberIds());
    }

    void SetupAnisotropyExtension(std::array<int, 2> fiber_ids)
    {
      anisotropyExtension_ =
          std::make_unique<MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension>(3, fiber_ids);

      anisotropy_.RegisterAnisotropyExtension(*anisotropyExtension_);

      anisotropy_.SetNumberOfGaussPoints(2);

      // Setup element fibers
      anisotropy_.SetGaussPointFibers(gpFibers_);
    }

    [[nodiscard]] int GetGaussPoint() const { return std::get<1>(GetParam()); }

    [[nodiscard]] std::array<int, 2> GetFiberIds() const { return std::get<0>(GetParam()); }

    MAT::Anisotropy anisotropy_;
    std::unique_ptr<MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension> anisotropyExtension_;

    std::vector<std::vector<CORE::LINALG::Matrix<3, 1>>> gpFibers_;
    std::vector<CORE::LINALG::Matrix<3, 3>> gpTensors_;
    std::vector<CORE::LINALG::Matrix<6, 1>> gpTensors_stress_;
    std::vector<double> gpScalarProducts_;
  };

  TEST_P(CoupAnisoExpoShearGaussPointFibersTest, GetScalarProduct)
  {
    EXPECT_NEAR(anisotropyExtension_->GetScalarProduct(GetGaussPoint()),
        gpScalarProducts_[GetGaussPoint()], 1e-10);
  }

  TEST_P(CoupAnisoExpoShearGaussPointFibersTest, GetStructuralTensor)
  {
    BACI_EXPECT_NEAR(anisotropyExtension_->GetStructuralTensor(GetGaussPoint()),
        gpTensors_[GetGaussPoint()], 1e-10);
  }

  TEST_P(CoupAnisoExpoShearGaussPointFibersTest, GetStructuralTensorStress)
  {
    BACI_EXPECT_NEAR(anisotropyExtension_->GetStructuralTensor_stress(GetGaussPoint()),
        gpTensors_stress_[GetGaussPoint()], 1e-10);
  }

  INSTANTIATE_TEST_SUITE_P(GaussPoints, CoupAnisoExpoShearGaussPointFibersTest,
      ::testing::Combine(::testing::Values(std::array<int, 2>({0, 1}), std::array<int, 2>({0, 2}),
                             std::array<int, 2>({1, 2})),
          ::testing::Values(0, 1)));
}  // namespace
